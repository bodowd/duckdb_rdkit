#include "duckdb/common/types.hpp"
#include "mol_descriptors.hpp"
#include "sdf_scanner/sdf_functions.hpp"

#define DUCKDB_EXTENSION_MAIN
#include "cast.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
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

static void LoadInternal(ExtensionLoader &loader) {
  duckdb_rdkit::RegisterTypes(loader);
  duckdb_rdkit::RegisterCasts(loader);
  duckdb_rdkit::RegisterFormatFunctions(loader);
  duckdb_rdkit::RegisterCompareFunctions(loader);
  duckdb_rdkit::RegisterDescriptorFunctions(loader);

  for (auto &fun : SDFFunctions::GetTableFunctions()) {
    loader.RegisterFunction(fun);
  }

  // SDF replacement scan
  auto &instance = loader.GetDatabaseInstance();
  auto &config = DBConfig::GetConfig(instance);
  config.replacement_scans.emplace_back(SDFFunctions::ReadSDFReplacement);
}

void DuckdbRdkitExtension::Load(ExtensionLoader &loader) {
  LoadInternal(loader);
}

std::string DuckdbRdkitExtension::Name() { return "duckdb_rdkit"; }

} // namespace duckdb

#ifdef DUCKDB_BUILD_LOADABLE_EXTENSION
extern "C" {

DUCKDB_CPP_EXTENSION_ENTRY(duckdb_rdkit, loader) {
  duckdb::LoadInternal(loader);
}

DUCKDB_EXTENSION_API const char *duckdb_rdkit_version() {
  return duckdb::DuckDB::LibraryVersion();
}
}
#endif

#ifndef DUCKDB_EXTENSION_MAIN
#error DUCKDB_EXTENSION_MAIN not defined
#endif
