#include "duckdb/common/types.hpp"
#include "mol_descriptors.hpp"
#include "sdf_scanner/sdf_functions.hpp"

#define DUCKDB_EXTENSION_MAIN
#include "cast.hpp"
#include "duckdb/main/extension_util.hpp"
#include "duckdb_rdkit_extension.hpp"
#include "mol_compare.hpp"
#include "mol_formats.hpp"
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

static void LoadInternal(DatabaseInstance &instance) {
  duckdb_rdkit::RegisterTypes(instance);
  duckdb_rdkit::RegisterCasts(instance);
  duckdb_rdkit::RegisterFormatFunctions(instance);
  duckdb_rdkit::RegisterCompareFunctions(instance);
  duckdb_rdkit::RegisterDescriptorFunctions(instance);

  for (auto &fun : SDFFunctions::GetTableFunctions()) {
    ExtensionUtil::RegisterFunction(instance, fun);
  }

  // SDF replacement scan
  auto &config = DBConfig::GetConfig(instance);
  config.replacement_scans.emplace_back(SDFFunctions::ReadSDFReplacement);
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
