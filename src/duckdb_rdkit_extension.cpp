#include "duckdb/common/types.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "mol_descriptors.hpp"
#include "sdf_scanner/sdf_functions.hpp"

#define DUCKDB_EXTENSION_MAIN
#include "cast.hpp"
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
  RegisterTypes(loader);
  RegisterCasts(loader);
  RegisterFormatFunctions(loader);
  RegisterCompareFunctions(loader);
  RegisterDescriptorFunctions(loader);

  for (auto &fun : SDFFunctions::GetTableFunctions()) {
    loader.RegisterFunction(fun);
  }

  // SDF replacement scan
  auto &config = DBConfig::GetConfig(loader.GetDatabaseInstance());
  config.replacement_scans.emplace_back(SDFFunctions::ReadSDFReplacement);
}

void DuckdbRdkitExtension::Load(ExtensionLoader &loader) {
  LoadInternal(loader);
}

std::string DuckdbRdkitExtension::Name() { return "duckdb_rdkit"; }

} // namespace duckdb

DUCKDB_CPP_EXTENSION_ENTRY(duckdb_rdkit, loader) {
  duckdb::LoadInternal(loader);
}
