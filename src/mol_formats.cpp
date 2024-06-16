#include "common.hpp"
#include "duckdb/common/string_util.hpp"
#include "types.hpp"
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace duckdb_rdkit {
// Expects a SMILES string and returns a RDKit pickled molecule
RDKit::ROMol *rdkit_mol_from_smiles(std::string s) {
  std::string smiles = s;
  RDKit::ROMol *mol;
  try {
    mol = RDKit::SmilesToMol(smiles);
  } catch (std::exception &e) {
    std::string msg = StringUtil::Format("%s", typeid(e).name());
    // not sure if this is the right way to throw an error in duckdb
    throw Exception(msg);
  }

  if (mol) {
    return mol;
  } else {
    std::string msg = StringUtil::Format("Could not convert %s to mol", smiles);
    throw Exception(msg);
  }
}

// Serialize a molecule to binary using RDKit's MolPickler
std::string rdkit_mol_to_binary_mol(const RDKit::ROMol mol) {
  std::string buf;
  try {
    RDKit::MolPickler::pickleMol(mol, buf);
  } catch (...) {
    std::string msg = "Could not serialize mol to binary";
    throw Exception(msg);
  }
  return buf;
}

// Deserialize a binary mol to RDKit mol
RDKit::ROMol rdkit_binary_mol_to_mol(std::string bmol) {
  RDKit::ROMol mol;
  RDKit::MolPickler::molFromPickle(bmol, mol);
  return mol;
}

// Expects an RDKit pickled molecule and returns the SMILES of the molecule
std::string rdkit_mol_to_smiles(RDKit::ROMol mol) {
  std::string smiles = RDKit::MolToSmiles(mol);
  return smiles;
}

// An extension function callable from duckdb
// converts a serialized RDKit molecule to SMILES
//
// If there is a table mols with a column of type Mol
//
//
// select mol_to_smiles(*) from mols;
//
void mol_to_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, string_t>(
      smiles, result, count, [&](string_t bmol) {
        auto mol = rdkit_binary_mol_to_mol(bmol.GetString());
        auto smiles = rdkit_mol_to_smiles(mol);
        return StringVector::AddString(result, smiles);
      });
}

// An extension function callable from duckdb
// converts a SMILES all the way to the serialized version of the RDKit mol
// returns NULL if conversion fails
//
// select mol_from_smiles('CC');
//
void mol_from_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles = args.data[0];
  auto count = args.size();

  UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
      smiles, result, count,
      [&](string_t smiles, ValidityMask &mask, idx_t idx) {
        try {
          auto mol = rdkit_mol_from_smiles(smiles.GetString());
          auto pickled_mol = rdkit_mol_to_binary_mol(*mol);
          return StringVector::AddString(result, pickled_mol);
        } catch (...) {
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
