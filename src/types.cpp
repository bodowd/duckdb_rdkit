#include "duckdb/common/types.hpp"
#include "duckdb/main/database.hpp"

namespace duckdb {

LogicalType Mol() {
  auto blob_type = LogicalType(LogicalTypeId::BLOB);
  blob_type.SetAlias("Mol");
  return blob_type;
}

void RegisterTypes(ExtensionLoader &loader) {
  loader.RegisterType("Mol", Mol());
}

} // namespace duckdb
