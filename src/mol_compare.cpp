#include "mol_compare.hpp"
#include "common.hpp"
#include "duckdb/common/enums/vector_type.hpp"
#include "duckdb/common/helper.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/unique_ptr.hpp"
#include "duckdb/execution/expression_executor_state.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <cstdio>
#include <memory>

namespace duckdb {

// credit: code is from chemicalite
// https://github.com/rvianello/chemicalite
// See mol_search.test for an example of
// a molecule which can return false negative, if the SMILES is different from
// the query
bool mol_cmp(std::string m1_bmol, std::string m2_bmol) {
  std::unique_ptr<RDKit::ROMol> m1(new RDKit::ROMol());
  std::unique_ptr<RDKit::ROMol> m2(new RDKit::ROMol());

  RDKit::MolPickler::molFromPickle(m1_bmol, *m1);
  RDKit::MolPickler::molFromPickle(m2_bmol, *m2);

  // credit: code is from chemicalite
  // https://github.com/rvianello/chemicalite
  // See mol_search.test for an example of
  // a molecule which can return false negative, if the SMILES is different
  // from the query if m1 is substruct of m2 and m2 is substruct of m1,
  // likely to be the same molecule
  RDKit::MatchVectType matchVect;
  bool recursion_possible = false;
  bool do_chiral_match = false; /* FIXME: make configurable getDoChiralSSS(); */
  bool ss1 = RDKit::SubstructMatch(*m1, *m2, matchVect, recursion_possible,
                                   do_chiral_match);
  bool ss2 = RDKit::SubstructMatch(*m2, *m1, matchVect, recursion_possible,
                                   do_chiral_match);
  if (ss1 && !ss2) {
    return false;
  } else if (!ss1 && ss2) {
    return false;
  }

  // the above can still fail in some chirality cases
  std::string smi1 = RDKit::MolToSmiles(*m1, do_chiral_match);
  std::string smi2 = RDKit::MolToSmiles(*m2, do_chiral_match);
  return smi1 == smi2;
}

static void is_exact_match(DataChunk &args, ExpressionState &state,
                           Vector &result) {
  D_ASSERT(args.ColumnCount() == 2);
  // args.data[i] is a FLAT_VECTOR
  auto &left = args.data[0];
  auto &right = args.data[1];

  BinaryExecutor::Execute<string_t, string_t, bool>(
      left, right, result, args.size(),
      [&](string_t &left_umbra_blob, string_t &right_umbra_blob) {
        auto left = umbra_mol_t(left_umbra_blob);
        auto right = umbra_mol_t(right_umbra_blob);

        // The prefix of a umbra_mol contains a bit vector for substructure
        // screens. We also use this to check exact match. If the molecules
        // being compared do not have the same substructures marked by the
        // dalke_fp, they cannot be an exact match
        if (memcmp(left.GetPrefix(), right.GetPrefix(),
                   umbra_mol_t::PREFIX_BYTES) != 0) {
          return false;
        };

        // otherwise, do the more extensive check with rdkit
        return mol_cmp(left.GetBinaryMol(), right.GetBinaryMol());
      });
}

bool _is_substruct(umbra_mol_t target, umbra_mol_t query,
                   RDKit::MatchVectType match_vect, bool recursion_possible,
                   bool do_chiral_match) {

  RDKit::ROMol left_mol;
  RDKit::ROMol right_mol;
  RDKit::MolPickler::molFromPickle(target.GetBinaryMol(), left_mol);
  RDKit::MolPickler::molFromPickle(target.GetBinaryMol(), right_mol);

  return RDKit::SubstructMatch(left_mol, right_mol, match_vect,
                               recursion_possible, do_chiral_match);
}

bool _is_substruct(umbra_mol_t target, RDKit::ROMol cached_mol,
                   RDKit::MatchVectType match_vect, bool recursion_possible,
                   bool do_chiral_match) {

  RDKit::ROMol left_mol;
  RDKit::MolPickler::molFromPickle(target.GetBinaryMol(), left_mol);

  return RDKit::SubstructMatch(left_mol, cached_mol, match_vect,
                               recursion_possible, do_chiral_match);
}

void IsSubstructLocalState::UpdateCache(string_t umbra_blob) {
  if (cached_mol == nullptr) {
    umbra_mol_t umbra_mol(umbra_blob);
    cached_mol = make_uniq<RDKit::ROMol>();
    RDKit::MolPickler::molFromPickle(umbra_mol.GetBinaryMol(), *cached_mol);
  }
}

unique_ptr<FunctionLocalState>
IsSubstructLocalState::Init(ExpressionState &state,
                            const BoundFunctionExpression &expr,
                            FunctionData *bind_data) {
  return make_uniq<IsSubstructLocalState>();
}

static void is_substruct(DataChunk &args, ExpressionState &state,
                         Vector &result) {
  D_ASSERT(args.ColumnCount() == 2);
  // args.data[i] is a FLAT_VECTOR
  auto &left = args.data[0];
  auto &right = args.data[1];

  auto &local_state = ExecuteFunctionState::GetFunctionState(state)
                          ->Cast<IsSubstructLocalState>();
  if (right.GetVectorType() == VectorType::CONSTANT_VECTOR) {
    auto right_data = ConstantVector::GetData<string_t>(right);
    local_state.UpdateCache(*right_data);
  }

  // copied from chemicalite
  RDKit::MatchVectType match_vect;
  bool recursion_possible = true;
  bool do_chiral_match = false; /* FIXME: make configurable getDoChiralSSS(); */

  BinaryExecutor::Execute<string_t, string_t, bool>(
      left, right, result, args.size(),
      [&](string_t &left_umbra_blob, string_t &right_umbra_blob) {
        auto query_umbra_mol = umbra_mol_t(right_umbra_blob);
        auto target_umbra_mol = umbra_mol_t(left_umbra_blob);
        // if the fragment exists in the query but not in the target,
        // there is no way for a match. This only works in one direction
        //
        // If the fragment exists in the target, but not the query, it is still
        // possible there is something in the query that matches the target, but
        // is not captured in the dalke fingerprint
        //
        // If all fragments that are on in the query are also on in the target,
        // this does not mean that the query is a substructure. It is possible
        // that there is something in the query not captured in the fingerprint
        // that is present in the query, but not in the target. For example,
        // if the query has NCCCCCCCC, and the target has the N bit set,
        // but it could be that the target is only NC
        //
        // It is only possible to short-circuit in the false case, not in the
        // true case
        auto q_prefix = query_umbra_mol.GetPrefixAsInt();
        auto t_prefix = target_umbra_mol.GetPrefixAsInt();
        // The 4 byte prefix in string_t is inlined. This is very fast to check.
        // If the first 4 bytes are a match, we need the rest of the dalke fp to
        // further check. This requires chasing a
        // pointer to the data of which the next 4 bytes of the dalke fp is at
        // the front of and is a little more expensive, so we only do this
        // if the first 4 bytes do not rule out the possiblity of being a
        // substructure match
        if ((q_prefix & t_prefix) == q_prefix) {
          auto q_dalke_fp = query_umbra_mol.GetDalkeFP();
          auto t_dalke_fp = target_umbra_mol.GetDalkeFP();
          // query might be substructure of the target -- run a substructure
          // match on the molecule objects
          if ((q_dalke_fp & t_dalke_fp) == q_dalke_fp) {
            if (local_state.cached_mol) {
              return _is_substruct(target_umbra_mol, *local_state.cached_mol,
                                   match_vect, recursion_possible,
                                   do_chiral_match);
            } else {
              // cached_mol could be empty if the right side (query molecule)
              // is not a constant, for example if the query compares two Mol
              // columns. In this case, we need to convert the query molecule
              // on every row because it possibly changes every row
              return _is_substruct(target_umbra_mol, query_umbra_mol,
                                   match_vect, recursion_possible,
                                   do_chiral_match);
            }
          }
        }
        return false;
      });
}

void RegisterCompareFunctions(ExtensionLoader &loader) {
  ScalarFunctionSet set("is_exact_match");
  // left type and right type
  set.AddFunction(
      ScalarFunction({Mol(), Mol()}, LogicalType::BOOLEAN, is_exact_match));
  loader.RegisterFunction(set);

  ScalarFunctionSet set_is_substruct("is_substruct");
  ScalarFunction is_substruct_func("is_substruct", {Mol(), Mol()},
                                   LogicalType::BOOLEAN, is_substruct);
  is_substruct_func.init_local_state = IsSubstructLocalState::Init;
  set_is_substruct.AddFunction(is_substruct_func);
  loader.RegisterFunction(set_is_substruct);
}

} // namespace duckdb
