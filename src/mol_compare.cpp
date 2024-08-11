#include "common.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/execution/expression_executor_state.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension_util.hpp"
#include "mol_formats.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <cstdio>
#include <memory>

namespace duckdb_rdkit {

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

        // The first 32 bits are the count prefix
        // The exact match only needs the count prefix
        if (memcmp(left.GetPrefix(), right.GetPrefix(),
                   umbra_mol_t::COUNT_PREFIX_BYTES) != 0) {
          return false;
        };

        // otherwise, do the more extensive check with rdkit
        return mol_cmp(left.GetBinaryMol(), right.GetBinaryMol());
      });
}

bool _is_substruct(umbra_mol_t target, umbra_mol_t query) {
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

  auto query_dalke_fp = query.GetDalkeFPBitset();
  auto target_dalke_fp = target.GetDalkeFPBitset();

  for (auto i = 0; i < query_dalke_fp.size(); i++) {
    if (query_dalke_fp[i] && !target_dalke_fp[i]) {
      return false;
    }
  }

  // otherwise, run a substructure match on the molecule objects
  std::unique_ptr<RDKit::ROMol> left_mol(new RDKit::ROMol());
  std::unique_ptr<RDKit::ROMol> right_mol(new RDKit::ROMol());

  RDKit::MolPickler::molFromPickle(target.GetBinaryMol(), *left_mol);
  RDKit::MolPickler::molFromPickle(query.GetBinaryMol(), *right_mol);

  // copied from chemicalite
  RDKit::MatchVectType matchVect;
  bool recursion_possible = true;
  bool do_chiral_match = false; /* FIXME: make configurable getDoChiralSSS(); */
  return RDKit::SubstructMatch(*left_mol, *right_mol, matchVect,
                               recursion_possible, do_chiral_match);
}

static void is_substruct(DataChunk &args, ExpressionState &state,
                         Vector &result) {
  D_ASSERT(args.ColumnCount() == 2);
  // args.data[i] is a FLAT_VECTOR
  auto &left = args.data[0];
  auto &right = args.data[1];

  BinaryExecutor::Execute<string_t, string_t, bool>(
      left, right, result, args.size(),
      [&](string_t &left_umbra_blob, string_t &right_umbra_blob) {
        auto left_umbra_mol = umbra_mol_t(left_umbra_blob);
        auto right_umbra_mol = umbra_mol_t(right_umbra_blob);

        return _is_substruct(left_umbra_mol, right_umbra_mol);
      });
}

void RegisterCompareFunctions(DatabaseInstance &instance) {
  ScalarFunctionSet set("is_exact_match");
  // left type and right type
  set.AddFunction(ScalarFunction({duckdb_rdkit::Mol(), duckdb_rdkit::Mol()},
                                 LogicalType::BOOLEAN, is_exact_match));
  ExtensionUtil::RegisterFunction(instance, set);

  ScalarFunctionSet set_is_substruct("is_substruct");
  set_is_substruct.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol(), duckdb_rdkit::Mol()},
                     LogicalType::BOOLEAN, is_substruct));
  ExtensionUtil::RegisterFunction(instance, set_is_substruct);
}

} // namespace duckdb_rdkit
