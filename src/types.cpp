// This file defines the new types in a way for duckdb to recognize

#include "types.hpp"
#include "common.hpp"
#include "duckdb/main/database.hpp"

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

std::vector<char> UmbraMol::serialize() {
  size_t total_size = HEADER_SIZE + sizeof(bmol_size) + bmol_size;
  std::vector<char> buffer;
  size_t offset = 0;

  // Copy each member to the buffer
  std::memcpy(buffer.data() + offset, &num_atoms, NUM_ATOMS_BYTES);
  offset += NUM_ATOMS_BYTES;
  std::memcpy(buffer.data() + offset, &num_bonds, NUM_BONDS_BYTES);
  offset += NUM_BONDS_BYTES;
  std::memcpy(buffer.data() + offset, &amw, AMW_BYTES);
  offset += AMW_BYTES;
  std::memcpy(buffer.data() + offset, &num_rings, NUM_RINGS_BYTES);
  offset += NUM_RINGS_BYTES;
  std::memcpy(buffer.data() + offset, &bmol_size, PTR_SIZE_BYTES);
  offset += PTR_SIZE_BYTES;
  std::memcpy(buffer.data() + offset, &bmol, bmol_size);

  return buffer;
}

void RegisterTypes(DatabaseInstance &instance) {
  // Register Mol type
  ExtensionUtil::RegisterType(instance, "Mol", duckdb_rdkit::Mol());
}

} // namespace duckdb_rdkit
