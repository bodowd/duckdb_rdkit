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
  static constexpr idx_t COUNT_PREFIX_BYTES = 4 * sizeof(char);
  static constexpr idx_t DALKE_FP_PREFIX_BYTES = 8 * sizeof(char);
  // 27 bits for the counts -- closest uint is 32 bits
  // 55 bits for the dalke_fp -- closest uint is 64 bits
  // using the uint32 and 64 makes concatenating stuff easier
  // as well, maybe this could be optimized and not waste
  // so much space
  // so 4 bytes + 8 byes = 12 bytes for prefix length
  static constexpr idx_t PREFIX_BYTES =
      COUNT_PREFIX_BYTES + DALKE_FP_PREFIX_BYTES;
  static constexpr idx_t INLINE_BYTES = 12 * sizeof(char);
  static constexpr idx_t MAX_STRING_SIZE = NumericLimits<uint32_t>::Maximum();
  static constexpr idx_t PREFIX_LENGTH = PREFIX_BYTES;

  umbra_mol_t() = default;

  // umbra_mol_t is a data type used for the duckdb_rdkit extension and it
  // is similar to the string_t type of duckdb.
  //
  // It has a prefix field which contains some calculated data from the
  // molecule, and a pointer to the binary molecule.
  // The prefix is used to short-circuit comparison operations; if there
  // is no chance for a match, it bails out and reduces work done by the system
  // by doing further more expensive work.
  //
  // In order to manage the data pointed to by the pointer when
  // it transitions between memory and the disk, pointer swizzling is required.
  // When duckdb sees a PhysicalType::VARCHAR, it can handle the pointer
  // swizzling.
  // The string_t is the memory representation of VARCHAR, and so the pointer
  // swizzling relies on string_t.
  //
  // umbra_mol_t is converted to and from string_t so that the rest of the
  // duckdb internals can handle the pointer swizzling, and whatever else
  // needs to happen.
  // string_t is the interface for umbra_mol_t to the rest of the duckdb system.
  // Or perhaps it can be thought of as the intermediate representation.
  umbra_mol_t(string_t buffer) {
    value.length = buffer.GetSize();
    memset(value.prefix, 0, PREFIX_LENGTH);
    memcpy(&value.prefix, buffer.GetData(), PREFIX_LENGTH);
    value.ptr = buffer.GetData();

    D_ASSERT(value.ptr == buffer.GetData());
  }

  std::bitset<64> GetDalkeFPBitset() {
    uint64_t int_fp = 0;
    // make sure to copy from value.prefix and not &value.prefix!
    // &value.prefix is not the data itself, but it is the address
    // of the prefix member of the struct.
    // Then adding 4 bytes to that gives an incorrect
    // memory address.
    // Instead, we want the data itself. value.prefix gives the starting
    // address of the char[] in value.prefix, which is the data itself.
    // Then we move the pointer 4 bytes, to go past the count prefix to the
    // start of the dalke_fp bits
    std::memcpy(&int_fp, value.prefix + COUNT_PREFIX_BYTES,
                DALKE_FP_PREFIX_BYTES);
    std::bitset<64> fp(int_fp);

    return fp;
  }

  const char *GetPrefix() { return value.prefix; }

  uint32_t GetBinaryMolSize() { return value.length - PREFIX_LENGTH; }

  std::string GetBinaryMol() {
    idx_t bmol_size = value.length - PREFIX_LENGTH;
    std::string buffer;
    buffer.resize(bmol_size);
    if (value.ptr && value.length > PREFIX_LENGTH) {
      memcpy(&buffer[0], &value.ptr[PREFIX_LENGTH], bmol_size);
    }
    return buffer;
  }

  idx_t GetSize() const { return value.length; }

  const char *GetData() const { return value.ptr; }

  std::string GetString() const { return std::string(GetData(), GetSize()); }

private:
  struct {
    uint32_t length;
    char prefix[PREFIX_LENGTH];
    const char *ptr;
  } value;
};

} // namespace duckdb_rdkit
