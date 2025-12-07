// This file defines the new types in a way for duckdb to recognize

#include "types.hpp"
#include "common.hpp"
#include "duckdb/common/constants.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/main/database.hpp"

namespace duckdb_rdkit {

LogicalType Mol() {
  auto blob_type = LogicalType(LogicalTypeId::BLOB);
  blob_type.SetAlias("Mol");
  return blob_type;
}

void RegisterTypes(ExtensionLoader &loader) {
  // Register Mol type
  loader.RegisterType("Mol", Mol());
}

} // namespace duckdb_rdkit
