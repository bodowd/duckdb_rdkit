#define DUCKDB_EXTENSION_MAIN

#include "duckdb_rdkit_extension.hpp"
#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension_util.hpp"
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>

namespace duckdb {

inline void DuckdbRdkitScalarFun(DataChunk &args, ExpressionState &state,
                                 Vector &result) {
  auto &name_vector = args.data[1];
  // try {
  //   std::unique_ptr<RDKit::ROMol> mol;
  //   mol.reset(RDKit::SmilesToMol(name_vector.ToString()));
  //   std::string smiles = RDKit::MolToSmiles(*mol);
  // } catch (...) {
  //   return;
  // };

  UnaryExecutor::Execute<string_t, string_t>(
      name_vector, result, args.size(), [&](string_t name) {
        return StringVector::AddString(result, "DuckdbRdkit " +
                                                   name.GetString() + " 🐥");
        ;
      });
}

static void LoadInternal(DatabaseInstance &instance) {
  // Register a scalar function
  auto duckdb_rdkit_scalar_function =
      ScalarFunction("duckdb_rdkit", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, DuckdbRdkitScalarFun);
  ExtensionUtil::RegisterFunction(instance, duckdb_rdkit_scalar_function);

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
