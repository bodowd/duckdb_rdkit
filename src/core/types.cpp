// This file defines the new types in a way for duckdb to recognize

#include "core/types.hpp"

#include "common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/operator/cast_operators.hpp"
#include "duckdb/common/unique_ptr.hpp"
#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "duckdb/parser/parsed_data/create_type_info.hpp"
#include <GraphMol/MolPickler.h>
#include <exception>
#include <string>

namespace duckdb_rdkit {
namespace core {

// Serialize a molecule to binary using RDKit's MolPickler
std::string mol_to_binary_mol(const RDKit::ROMol &mol, int *rc) {
  std::string buf;
  try {
    RDKit::MolPickler::pickleMol(mol, buf, RDKit::PicklerOps::AllProps);
  } catch (...) {

    std::string msg = "Could not serialize mol to binary";
    // NOTE: not sure if this is the right way to throw an error in duckdb
    HandleCastError::AssignError(msg, &msg);
  }
  return buf;
}

// The mol object holds a lot of information about the atoms and bonds, etc.
// of the molecule. This information is used for further operations on the
// molecule. The Mol object will be serialized to binary using RDKit's
// MolPickler and stored as binary. Binary is stored as BLOB in duckdb
LogicalType RDKitTypes::Mol() {
  auto type = LogicalType::BLOB;
  return type;
}

void RDKitTypes::Register(DatabaseInstance &db) {
  // Mol
  ExtensionUtil::RegisterType(db, "Mol", RDKitTypes::Mol());
}

} // namespace core
} // namespace duckdb_rdkit
