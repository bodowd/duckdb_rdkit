#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/execution/expression_executor_state.hpp"
#include "duckdb/function/function_set.hpp"
#include <exception>
#define DUCKDB_EXTENSION_MAIN

#include "duckdb/common/exception.hpp"
#include "duckdb/common/operator/cast_operators.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension_util.hpp"
#include "duckdb_rdkit_extension.hpp"
#include "types.hpp"
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>
#include <string>

namespace duckdb {

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

// Expects a SMILES string and returns a RDKit pickled molecule
static RDKit::ROMol *rdkit_mol_from_smiles(std::string s) {
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

// Expects an RDKit pickled molecule and returns the SMILES of the molecule
static std::string rdkit_mol_to_smiles(RDKit::ROMol mol) {
  std::string smiles = RDKit::MolToSmiles(mol);
  return smiles;
}

// An extension function callable from duckdb
// converts a SMILES all the way to the serialized version of the RDKit mol
//
// select mol_from_smiles('CC');
//
void mol_from_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, string_t>(
      smiles, result, count, [&](string_t smiles) {
        auto mol = rdkit_mol_from_smiles(smiles.GetString());
        auto pickled_mol = rdkit_mol_to_binary_mol(*mol);
        return StringVector::AddString(result, pickled_mol);
      });
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

// This enables the user to insert into a Mol column by just writing the SMILES
// Duckdb will try to convert the string to a rdkit mol
// This is consistent with the RDKit Postgres cartridge behavior
void VarcharToMol(Vector &source, Vector &result, idx_t count) {
  UnaryExecutor::Execute<string_t, string_t>(
      source, result, count, [&](string_t smiles) {
        auto mol = rdkit_mol_from_smiles(smiles.GetString());
        auto pickled_mol = rdkit_mol_to_binary_mol(*mol);
        return StringVector::AddString(result, pickled_mol);
      });
}

bool VarcharToMolCast(Vector &source, Vector &result, idx_t count,
                      CastParameters &parameters) {
  VarcharToMol(source, result, count);
  return true;
}

void MolToVarchar(Vector &source, Vector &result, idx_t count) {
  UnaryExecutor::Execute<string_t, string_t>(
      source, result, count, [&](string_t bmol) {
        auto mol = rdkit_binary_mol_to_mol(bmol.GetString());
        auto smiles = rdkit_mol_to_smiles(mol);
        return StringVector::AddString(result, smiles);
      });
}

bool MolToVarcharCast(Vector &source, Vector &result, idx_t count,
                      CastParameters &parameters) {
  MolToVarchar(source, result, count);
  return true;
}

static void LoadInternal(DatabaseInstance &instance) {
  // Register Mol type
  ExtensionUtil::RegisterType(instance, "Mol", duckdb_rdkit::Mol());

  // Register scalar functions

  ScalarFunctionSet mol_from_smiles_set("mol_from_smiles");
  mol_from_smiles_set.AddFunction(ScalarFunction(
      {LogicalType::VARCHAR}, duckdb_rdkit::Mol(), mol_from_smiles));
  ExtensionUtil::RegisterFunction(instance, mol_from_smiles_set);

  ScalarFunctionSet mol_to_smiles_set("mol_to_smiles");
  mol_to_smiles_set.AddFunction(ScalarFunction(
      {duckdb_rdkit::Mol()}, LogicalType::VARCHAR, mol_to_smiles));
  ExtensionUtil::RegisterFunction(instance, mol_to_smiles_set);

  ExtensionUtil::RegisterCastFunction(instance, LogicalType::VARCHAR,
                                      duckdb_rdkit::Mol(),
                                      BoundCastInfo(VarcharToMolCast), 1);

  ExtensionUtil::RegisterCastFunction(instance, duckdb_rdkit::Mol(),
                                      LogicalType::VARCHAR,
                                      BoundCastInfo(MolToVarcharCast), 1);
}

void DuckdbRdkitExtension::Load(DuckDB &db) { LoadInternal(*db.instance); }

std::string DuckdbRdkitExtension::Name() { return "duckdb_rdkit"; }

} // namespace duckdb

extern "C" {

DUCKDB_EXTENSION_API void duckdb_rdkit_init(duckdb::DatabaseInstance &db) {
  duckdb::DuckDB db_wrapper(db);
  db_wrapper.LoadExtension<duckdb::DuckdbRdkitExtension>();
}

DUCKDB_EXTENSION_API const char *duckdb_rdkit_version() {
  return duckdb::DuckDB::LibraryVersion();
}
}

#ifndef DUCKDB_EXTENSION_MAIN
#error DUCKDB_EXTENSION_MAIN not defined
#endif
