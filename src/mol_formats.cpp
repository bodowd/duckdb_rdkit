#include "mol_formats.hpp"
#include "common.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function_set.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace duckdb_rdkit {
// Expects a SMILES string and returns a RDKit pickled molecule
std::unique_ptr<RDKit::ROMol> rdkit_mol_from_smiles(std::string s) {
  std::string smiles = s;
  std::unique_ptr<RDKit::ROMol> mol;
  try {
    mol.reset(RDKit::SmilesToMol(smiles));
  } catch (std::exception &e) {
    std::string msg = StringUtil::Format("%s", typeid(e).name());
    throw InvalidInputException(msg);
  }

  if (mol) {
    return mol;
  } else {
    std::string msg = StringUtil::Format("Could not convert %s to mol", smiles);
    throw InvalidInputException(msg);
  }
}

// Serialize a molecule to binary using RDKit's MolPickler
std::string rdkit_mol_to_binary_mol(const RDKit::ROMol mol) {
  std::string buf;
  try {
    RDKit::MolPickler::pickleMol(mol, buf);
  } catch (...) {
    std::string msg = "Could not serialize mol to binary";
    throw InvalidInputException(msg);
  }
  return buf;
}

// Deserialize a binary mol to RDKit mol
std::unique_ptr<RDKit::ROMol> rdkit_binary_mol_to_mol(std::string bmol) {
  std::unique_ptr<RDKit::ROMol> mol(new RDKit::ROMol());
  RDKit::MolPickler::molFromPickle(bmol, *mol);

  return mol;
}

std::string rdkit_mol_to_smiles(RDKit::ROMol mol) {
  std::string smiles = RDKit::MolToSmiles(mol);
  return smiles;
}

void mol_to_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &bmol = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, string_t>(
      bmol, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto bmol = umbra_mol.GetBinaryMol();
        auto mol = rdkit_binary_mol_to_mol(bmol);
        auto smiles = rdkit_mol_to_smiles(*mol);
        return StringVector::AddString(result, smiles);
      });
}

void mol_from_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles = args.data[0];
  auto count = args.size();

  UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
      smiles, result, count,
      [&](string_t smiles, ValidityMask &mask, idx_t idx) {
        try {
          auto mol = rdkit_mol_from_smiles(smiles.GetString());

          auto res = get_umbra_mol_string(*mol);

          // IMPORTANT! StringVector::AddString needs to take a std::string
          // Using string_t::GetString() seems to mangle the data
          return StringVector::AddStringOrBlob(result, res);
        } catch (...) {
          std::cout << "WARNING: could not create molecule from SMILES\n"
                    << smiles.GetData() << std::endl;
          // printf("WARNING: could not create molecule from SMILES %s\n",
          //        smiles.GetData());
          mask.SetInvalid(idx);
          return string_t();
        }
      });
}

void RegisterFormatFunctions(DatabaseInstance &instance) {
  // Register scalar functions
  ScalarFunctionSet mol_from_smiles_set("mol_from_smiles");
  mol_from_smiles_set.AddFunction(ScalarFunction(
      {LogicalType::VARCHAR}, duckdb_rdkit::Mol(), mol_from_smiles));
  ExtensionUtil::RegisterFunction(instance, mol_from_smiles_set);

  ScalarFunctionSet mol_to_smiles_set("mol_to_smiles");
  mol_to_smiles_set.AddFunction(ScalarFunction(
      {duckdb_rdkit::Mol()}, LogicalType::VARCHAR, mol_to_smiles));
  ExtensionUtil::RegisterFunction(instance, mol_to_smiles_set);
}

} // namespace duckdb_rdkit
