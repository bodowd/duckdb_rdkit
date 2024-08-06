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
bool mol_cmp(const RDKit::ROMol &m1, const RDKit::ROMol &m2) {
  int res = m1.getNumAtoms() - m2.getNumAtoms();

  if (res) {
    return false;
  }

  res = m1.getNumBonds() - m2.getNumBonds();
  if (res) {
    return false;
  }

  res = int(RDKit::Descriptors::calcAMW(m1, false) -
            RDKit::Descriptors::calcAMW(m2, false) + .5);
  if (res) {
    return false;
  }

  res = m1.getRingInfo()->numRings() - m2.getRingInfo()->numRings();
  if (res) {
    return false;
  }

  // if m1 is substruct of m2 and m2 is substruct of m1, likely to be the same
  // molecule
  RDKit::MatchVectType matchVect;
  bool recursion_possible = false;
  bool do_chiral_match = false; /* FIXME: make configurable getDoChiralSSS(); */
  bool ss1 = RDKit::SubstructMatch(m1, m2, matchVect, recursion_possible,
                                   do_chiral_match);
  bool ss2 = RDKit::SubstructMatch(m2, m1, matchVect, recursion_possible,
                                   do_chiral_match);
  if (ss1 && !ss2) {
    return false;
  } else if (!ss1 && ss2) {
    return false;
  }

  // the above can still fail in some chirality cases
  std::string smi1 = RDKit::MolToSmiles(m1, do_chiral_match);
  std::string smi2 = RDKit::MolToSmiles(m2, do_chiral_match);
  return smi1 == smi2;
}

// runs exact match on two duckdb Flat vectors
static void is_exact_match(DataChunk &args, ExpressionState &state,
                           Vector &result) {
  D_ASSERT(args.ColumnCount() == 2);
  // args.data[i] is a FLAT_VECTOR
  auto &left = args.data[0];
  auto &right = args.data[1];

  BinaryExecutor::Execute<string_t, string_t, bool>(
      left, right, result, args.size(),
      [&](string_t &left_blob, string_t &right_blob) {
        // need to first deserialize the blob which is stored in the DB
        // NOTE: is it possible to run some of the search on the binary mol?
        // for example, comparing MW. Is this stored in the binary,
        // and is it possible to jump to the offset where that info is stored
        // and run the comparison of there?
        std::unique_ptr<RDKit::ROMol> left_mol(new RDKit::ROMol());
        std::unique_ptr<RDKit::ROMol> right_mol(new RDKit::ROMol());

        RDKit::MolPickler::molFromPickle(left_blob.GetString(), *left_mol);
        RDKit::MolPickler::molFromPickle(right_blob.GetString(), *right_mol);
        auto compare_result = mol_cmp(*left_mol, *right_mol);
        return compare_result;
      });
}

// bool umbra_mol_cmp(std::string m1_bmol, std::string m2_bmol) {
//
//   // otherwise, run a full check on the molecule objects
//   std::unique_ptr<RDKit::ROMol> left_mol(new RDKit::ROMol());
//   std::unique_ptr<RDKit::ROMol> right_mol(new RDKit::ROMol());
//
//   RDKit::MolPickler::molFromPickle(m1_bmol, *left_mol);
//   RDKit::MolPickler::molFromPickle(m2_bmol, *right_mol);
//
//   // experiment: log when the above check does not short circuit
//   // {
//   //   std::ofstream log_file("log_file.txt",
//   //                          std::ios_base::out | std::ios_base::app);
//   //   log_file << "left_mol: " << rdkit_mol_to_smiles(*left_mol) << ","
//   //            << "right_mol: " << rdkit_mol_to_smiles(*right_mol) <<
//   //            std::endl;
//   // }
//   return mol_cmp(*left_mol, *right_mol);
// }
bool umbra_mol_cmp(std::string m1_bmol, std::string m2_bmol) {
  // experiment: log when the above check does not short circuit
  // {
  //   std::ofstream log_file("log_file.txt",
  //                          std::ios_base::out | std::ios_base::app);
  //   log_file << "left_mol: " << rdkit_mol_to_smiles(*left_mol) << ","
  //            << "right_mol: " << rdkit_mol_to_smiles(*right_mol) <<
  //            std::endl;
  // }

  // otherwise, run a full check on the molecule objects
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

static void umbra_is_exact_match(DataChunk &args, ExpressionState &state,
                                 Vector &result) {
  D_ASSERT(args.ColumnCount() == 2);
  // args.data[i] is a FLAT_VECTOR
  auto &left = args.data[0];
  auto &right = args.data[1];

  BinaryExecutor::Execute<string_t, string_t, bool>(
      left, right, result, args.size(),
      [&](string_t &left_umbra_blob, string_t &right_umbra_blob) {
        std::cout << "umbra_is_exact_match: " << std::endl;
        std::cout << "\nleft_umbra_blob: " << std::endl;
        auto left_ptr = left_umbra_blob.GetData();
        for (auto i = 0; i < left_umbra_blob.GetSize(); i++) {
          printf("%02x ", static_cast<unsigned char>(left_ptr[i]));
        }
        std::cout << "\nright_umbra_blob: " << std::endl;
        auto right_ptr = right_umbra_blob.GetData();
        for (auto i = 0; i < right_umbra_blob.GetSize(); i++) {
          printf("%02x ", static_cast<unsigned char>(right_ptr[i]));
        }
        auto left = umbra_mol_t(left_umbra_blob);
        auto right = umbra_mol_t(right_umbra_blob);

        for (char b : left.GetString()) {
          printf("%02x ", static_cast<unsigned char>(b));
        }

        // first 27 bits of the prefix are the counts
        // That fits into a uint32_t
        // First copy the first 4 bytes of the prefix
        // Then get the 27 bits of the prefix
        uint32_t a_count_prefix;
        memcpy(&a_count_prefix, left.GetPrefix(),
               umbra_mol_t::COUNT_PREFIX_BYTES);

        uint32_t b_count_prefix;
        memcpy(&b_count_prefix, right.GetPrefix(),
               umbra_mol_t::COUNT_PREFIX_BYTES);
        // shift to the right to get the highest 27 bits
        // The counts prefix are packed all the way to the highest
        // bit, number 31 (counting from 0), so the lowest 5 bits are not
        // part of the count prefix
        // uint32_t a_count_prefix = a_prefix_32_bits >> 5;
        // uint32_t b_count_prefix = b_prefix_32_bits >> 5;

        // auto a_prefix = Load<uint32_t>(const_data_ptr_cast(a.GetPrefix()));
        // uint16_t b_prefix =
        // Load<uint32_t>(const_data_ptr_cast(b.GetPrefix()));

        if (a_count_prefix != b_count_prefix) {
          return false;
        }

        // otherwise, do the more extensive check with rdkit

        // std::string a_bmol;
        // std::string b_bmol;
        // memcpy(&a_bmol[0], left.GetBinaryMol(), left.GetBina());
        // memcpy(&b_bmol[0], right.GetBinaryMol(), right.GetBmolSize());
        //
        // return umbra_mol_cmp(a_bmol, b_bmol);

        return true;
      });
}

void RegisterCompareFunctions(DatabaseInstance &instance) {
  ScalarFunctionSet set("is_exact_match");
  // left type and right type
  set.AddFunction(ScalarFunction({duckdb_rdkit::Mol(), duckdb_rdkit::Mol()},
                                 LogicalType::BOOLEAN, is_exact_match));
  ExtensionUtil::RegisterFunction(instance, set);

  ScalarFunctionSet set_umbra_exact_match("umbra_is_exact_match");
  set_umbra_exact_match.AddFunction(
      ScalarFunction({duckdb_rdkit::UmbraMol(), duckdb_rdkit::UmbraMol()},
                     LogicalType::BOOLEAN, umbra_is_exact_match));
  ExtensionUtil::RegisterFunction(instance, set_umbra_exact_match);
}

} // namespace duckdb_rdkit
