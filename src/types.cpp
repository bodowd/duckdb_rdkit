// This file defines the new types in a way for duckdb to recognize

#include "types.hpp"
#include "common.hpp"
#include "duckdb/common/constants.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/extension_util.hpp"

namespace duckdb_rdkit {

LogicalType Mol() {
  auto blob_type = LogicalType(LogicalTypeId::BLOB);
  blob_type.SetAlias("Mol");
  return blob_type;
}

void RegisterTypes(DatabaseInstance &instance) {
  // Register Mol type
  ExtensionUtil::RegisterType(instance, "Mol", duckdb_rdkit::Mol());
}

} // namespace duckdb_rdkit
