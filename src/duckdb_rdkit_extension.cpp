#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/string_type.hpp"
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
//
std::string rdkit_mol_to_binary_mol(const RDKit::ROMol &mol) {
  std::string buf;
  try {
    RDKit::MolPickler::pickleMol(mol, buf, RDKit::PicklerOps::AllProps);
  } catch (...) {
    std::string msg = "Could not serialize mol to binary";
    // not sure if this is the right way to throw an error in duckdb
    // HandleCastError::AssignError(msg, &msg);
    throw Exception(msg);
  }
  return buf;
}

static std::string rdkit_mol_from_smiles(std::string s) {
  /* build the molecule blob repr from a text string */
  std::string smiles = s;
  std::unique_ptr<RDKit::ROMol> mol;

  try {
    mol.reset(RDKit::SmilesToMol(smiles));
  } catch (std::exception &e) {
    // std::string msg = StringUtil::Format("Could not convert %s to mol",
    // smiles);
    std::string msg = StringUtil::Format("%s", typeid(e).name());
    // not sure if this is the right way to throw an error in duckdb
    throw Exception(msg);
  }

  if (mol) {
    // serialize the mol
    std::string bmol = rdkit_mol_to_binary_mol(*mol);
    return bmol;
  } else {
    std::string msg = StringUtil::Format("Could not convert %s to mol", smiles);
    throw Exception(msg);
  }
}

inline void DuckdbRdkitScalarFun(DataChunk &args, ExpressionState &state,
                                 Vector &result) {
  auto &smiles_vector = args.data[0];
  UnaryExecutor::Execute<string_t, string_t>(
      smiles_vector, result, args.size(), [&](string_t s) {
        auto b = rdkit_mol_from_smiles(s.GetString());
        std::string msg = StringUtil::Format("%s", b);
        HandleCastError::AssignError(msg, &msg);
        return StringVector::AddString(result, "DuckdbRdkit " + msg + "🐥");
        ;
      });
}

void mol_from_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, string_t>(
      smiles, result, count, [&](string_t smiles) {
        auto binary_mol = rdkit_mol_from_smiles(smiles.GetString());
        return binary_mol;
      });
}

// void mol_to_smiles() {
//   std::unique_ptr<RDKit::ROMol> mol(arg_to_romol(arg, &rc));
// }

static void LoadInternal(DatabaseInstance &instance) {
  // Register a scalar function
  auto duckdb_rdkit_scalar_function =
      ScalarFunction("duckdb_rdkit", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, DuckdbRdkitScalarFun);
  ExtensionUtil::RegisterFunction(instance, duckdb_rdkit_scalar_function);

  ScalarFunctionSet mol_from_smiles_set("mol_from_smiles");
  mol_from_smiles_set.AddFunction(ScalarFunction(
      {LogicalType::VARCHAR}, duckdb_rdkit::Mol(), mol_from_smiles));
  ExtensionUtil::RegisterFunction(instance, mol_from_smiles_set);
  // ScalarFunctionSet as_wkb_function_set("ST_AsWKB");
  //
  // as_wkb_function_set.AddFunction(ScalarFunction(
  //     {GeoTypes::GEOMETRY()}, GeoTypes::WKB_BLOB(), GeometryAsWBKFunction,
  //     nullptr, nullptr, nullptr, GeometryFunctionLocalState::Init));
  //
  // ExtensionUtil::RegisterFunction(db, as_wkb_function_set);

  // Register another scalar function
  // auto duckdb_rdkit_openssl_version_scalar_function =
  //     ScalarFunction("duckdb_rdkit_openssl_version", {LogicalType::VARCHAR},
  //                    LogicalType::VARCHAR,
  //                    DuckdbRdkitOpenSSLVersionScalarFun);
  // ExtensionUtil::RegisterFunction(instance,
  //                                 duckdb_rdkit_openssl_version_scalar_function);
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
