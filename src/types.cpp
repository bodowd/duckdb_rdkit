// This file defines the new types in a way for duckdb to recognize

#include "types.hpp"
#include "common.hpp"
#include "duckdb/common/constants.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/extension_util.hpp"

namespace duckdb_rdkit {

// The mol object holds a lot of information about the atoms and bonds, etc.
// of the molecule. This information is used for further operations on the
// molecule. The Mol object will be serialized to binary using RDKit's
// MolPickler and stored as binary. Binary is stored as BLOB in duckdb
LogicalType Mol() {
  auto blob_type = LogicalType(LogicalTypeId::BLOB);
  blob_type.SetAlias("Mol");
  return blob_type;
}

LogicalType UmbraMol() {
  auto type = LogicalType(LogicalTypeId::BLOB);
  type.SetAlias("UmbraMol");
  return type;
}

void RegisterTypes(DatabaseInstance &instance) {
  // Register Mol type
  ExtensionUtil::RegisterType(instance, "Mol", duckdb_rdkit::Mol());
  ExtensionUtil::RegisterType(instance, "UmbraMol", duckdb_rdkit::UmbraMol());
}

} // namespace duckdb_rdkit
