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

// This is to generate the prefix and concatenate with the binary RDKit molecule
// so that it can then be sent to a string_t. Return the std::string because
// later the StringVector::AddStringOrBlob function takes a std::string, not
// string_t
std::string get_umbra_mol_string(const RDKit::ROMol &mol);

struct umbra_mol_t {
  static constexpr idx_t COUNT_PREFIX_BYTES = 4 * sizeof(char);
  static constexpr idx_t DALKE_FP_PREFIX_BYTES = 8 * sizeof(char);
  // 27 bits for the counts -- closest uint is 32 bits
  // 55 bits for the dalke_fp -- closest uint is 64 bits
  // using the uint32 and 64 makes concatenating stuff easier
  // for now, maybe this could be optimized and not waste
  // so much space
  // so 4 bytes + 8 byes = 12 bytes for prefix length
  static constexpr idx_t PREFIX_BYTES =
      COUNT_PREFIX_BYTES + DALKE_FP_PREFIX_BYTES;
  static constexpr idx_t INLINE_BYTES = 12 * sizeof(char);
  static constexpr idx_t HEADER_SIZE = sizeof(uint32_t) + PREFIX_BYTES;
  static constexpr idx_t MAX_STRING_SIZE = NumericLimits<uint32_t>::Maximum();
  static constexpr idx_t PREFIX_LENGTH = PREFIX_BYTES;
  static constexpr idx_t INLINE_LENGTH = INLINE_BYTES;

  umbra_mol_t() = default;

  // umbra_mol_t is a data type used for the duckdb_rdkit extension and it
  // is similar to the string_t type of duckdb.
  // umbra_mol_t keeps the string in a pointer and only follows the pointer
  // to the string when the fast checks using the prefix does not answer
  // rule things out in the exact and substructure match functions.
  // In order to manage the data pointed to by the pointer when it transitions
  // between memory and the disk, pointer swizzling is employed by duckdb, and
  // that implementation checks for a Physical Type VARCHAR. The string_t
  // is the memory representation of VARCHAR, and so the pointer swizzling
  // relies on string_t. Since this is internal to duckdb, we use string_t
  // as an intermediate between memory and disk for umbra_mol_t.
  //
  // when a varchar is to be converted to umbra_mol_t, we simply need to
  // generate the binary data string that contains the desired data for
  // UmbraMol, and then pass that to a string_t function. Duckdb handles the
  // rest there for getting it to disk and pointer swizzling.
  //
  // When we want to use that binary data, we also get it out of the system as
  // a string_t. Again, duckdb handles the pointer swizzling. Then we convert
  // string_t to an umbra_mol by extracting the bytes desired for the prefix,
  // and then setting the ptr field of umbra_mol_t to the ptr field of string_t
  //
  // string_t is the interface for umbra_mol_t to the rest of the duckdb system.
  // Or perhaps it can be thought of as the intermediate representation.
  umbra_mol_t(string_t buffer) {
    value.length = buffer.GetSize();
    memset(value.prefix, 0, PREFIX_LENGTH);
    memcpy(&value.prefix, buffer.GetPrefix(), PREFIX_LENGTH);
    // std::cout << "\nCONSTRUCTOR:" << std::endl;
    // for (char b : buffer.GetString()) {
    //   printf("%02x ", static_cast<unsigned char>(b));
    // }
    // std::cout << "\nBuffer size: " << buffer.GetSize() << std::endl;
    // std::cout << "Print from buffer pointer: " << std::endl;
    // auto p = buffer.GetData();
    // for (auto i = 0; i < buffer.GetSize(); i++) {
    //   printf("%02x ", static_cast<unsigned char>(p[i]));
    // }

    // std::cout << "Print from value.ptr: " << std::endl;
    value.ptr = buffer.GetData();
    // for (auto i = 0; i < buffer.GetSize(); i++) {
    //   printf("%02x ", static_cast<unsigned char>(value.ptr[i]));
    // }

    D_ASSERT(value.ptr == buffer.GetData());
    // std::cout << "value.ptr == buffer.GetData(): "
    //           << (value.ptr == buffer.GetData()) << std::endl;
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
    // this is the binary molecule size
    // the prefix will be calculated in the umbra_mol_t struct constructor
    uint32_t length;
    char prefix[PREFIX_LENGTH];
    const char *ptr;
  } value;
};

} // namespace duckdb_rdkit
