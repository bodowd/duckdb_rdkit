#pragma once
#include "common.hpp"
#include "duckdb/common/shared_ptr.hpp"
#include "duckdb/common/typedefs.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <cstdint>
#include <cstring>

namespace duckdb_rdkit {

// This is to generate the prefix and concatenate it with the binary RDKit
// molecule so that it can then be sent to a string_t. Return the std::string
// because later the StringVector::AddStringOrBlob function takes a std::string,
// not string_t.
std::string get_umbra_mol_string(const RDKit::ROMol &mol);

struct umbra_mol_t {
  // Use composition to add methods to the string_t
  // the umbra_mol is just a string_t under the hood.
  // There are only additional methods on top to extract the binary
  // Mol or get a certain prefix, for example.
  // This should not require a copy of the string_t
  string_t &string_t_umbra_mol;

  umbra_mol_t(string_t &buffer) : string_t_umbra_mol(buffer) {}

  // 55 bits for the dalke_fp -- closest uint is 64 bits
  // This is the entire Dalke FP size
  // 4 bytes are inlined by duckdb as the PREFIX in the underlying string_t
  // class And the remaining 4 bytes are in the beginning of the "string"
  // pointed to by the pointer in string_t
  static constexpr idx_t DALKE_FP_PREFIX_BYTES = 8 * sizeof(char);
  static constexpr idx_t MAX_STRING_SIZE = NumericLimits<uint32_t>::Maximum();
  static constexpr idx_t PREFIX_BYTES = string_t::PREFIX_BYTES;

  // umbra_mol_t is a data type used for the duckdb_rdkit extension and it
  // is a string_t type under the hood.
  //
  // It has a prefix field which contains some calculated data from the
  // molecule, and a pointer to the binary molecule.
  // The prefix is used to short-circuit comparison operations; if there
  // is no chance for a match, it bails out and reduces work done by the
  //  system
  // by doing further more expensive work.
  //
  // In order to manage the data pointed to by the pointer when
  // it transitions between memory and the disk, pointer swizzling is
  //     required.
  // When duckdb sees a PhysicalType::VARCHAR, it can handle the pointer
  // swizzling.
  // The string_t is the memory representation of VARCHAR, and so the
  // pointer swizzling relies on string_t.
  //
  // umbra_mol_t is converted to and from string_t so that the rest of the
  // duckdb internals can handle the pointer swizzling, and whatever else
  // needs to happen.
  // string_t is the interface for umbra_mol_t to the rest of the duckdb
  // system.
  // Or perhaps it can be thought of as the intermediate representation.

  // umbra_mol_t(string_t buffer) {
  //   value.binary_umbra_mol = buffer;
  //   value.length = buffer.GetSize();
  //   memset(value.prefix, 0, PREFIX_LENGTH);
  //   memcpy(&value.prefix, buffer.GetData(), PREFIX_LENGTH);
  //   value.ptr = buffer.GetData();
  //
  //   D_ASSERT(value.ptr == buffer.GetData());
  // }

  uint64_t GetDalkeFP() {
    uint64_t int_fp = 0;
    std::memcpy(&int_fp, string_t_umbra_mol.GetData(), DALKE_FP_PREFIX_BYTES);
    return int_fp;
  }

  // Return the prefix as a 4 byte int
  // Converts the underlying string_t prefix to 4 byte int to make it
  // easy to do bitwise operation
  uint32_t GetPrefixAsInt() {
    return Load<uint32_t>(const_data_ptr_cast(string_t_umbra_mol.GetPrefix()));
  }

  const char *GetPrefix() { return string_t_umbra_mol.GetPrefix(); }

  uint32_t GetBinaryMolSize() {
    return string_t_umbra_mol.GetSize() - DALKE_FP_PREFIX_BYTES;
  }

  std::string GetBinaryMol() {
    idx_t bmol_size = string_t_umbra_mol.GetSize() - DALKE_FP_PREFIX_BYTES;
    std::string buffer;
    buffer.resize(bmol_size);
    if (string_t_umbra_mol.GetData() &&
        string_t_umbra_mol.GetSize() > DALKE_FP_PREFIX_BYTES) {
      memcpy(&buffer[0], &string_t_umbra_mol.GetData()[DALKE_FP_PREFIX_BYTES],
             bmol_size);
    }
    return buffer;
  }

  idx_t GetSize() const { return string_t_umbra_mol.GetSize(); }

  const char *GetData() const { return string_t_umbra_mol.GetData(); }

  std::string GetString() const { return std::string(GetData(), GetSize()); }
};

} // namespace duckdb_rdkit
