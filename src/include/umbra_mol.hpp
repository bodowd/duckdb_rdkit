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

std::string count_prefix(uint32_t num_atoms, uint32_t num_bonds, uint32_t amw,
                         uint32_t num_rings, const std::string &binary_mol);

struct umbra_mol_t {
  static constexpr idx_t COUNT_PREFIX_BYTES = 4 * sizeof(char);
  static constexpr idx_t DALKE_FP_PREFIX_BYTES = 8 * sizeof(char);
  // 27 bits for the counts -- closest uint is 32 bits
  // 55 bits for the dalke_fp -- closest uint is 64 bits
  // using the uint32 and 64 makes concatenating stuff easier
  // for now, maybe this could be optimized and not waste
  // so much space
  // so 4 bytes + 8 byes = 12 bytes for prefix length
  static constexpr idx_t PREFIX_BYTES = 12 * sizeof(char);
  static constexpr idx_t INLINE_BYTES = 12 * sizeof(char);
  static constexpr idx_t HEADER_SIZE = sizeof(uint32_t) + PREFIX_BYTES;
  static constexpr idx_t MAX_STRING_SIZE = NumericLimits<uint32_t>::Maximum();
  static constexpr idx_t PREFIX_LENGTH = PREFIX_BYTES;
  static constexpr idx_t INLINE_LENGTH = INLINE_BYTES;

  umbra_mol_t() = default;

  umbra_mol_t(string_t buffer) {
    value.length = buffer.GetSize();
    memset(value.prefix, 0, PREFIX_LENGTH);
    memcpy(&value.prefix, buffer.GetPrefix(), PREFIX_LENGTH);
    // the lambda expression `[](char*){}` is a custom deleter for cleaning up
    // when the last owner of the shared_ptr is destroyed
    std::cout << "HERE: " << std::endl;
    for (char b : buffer.GetString()) {
      printf("%02x ", static_cast<unsigned char>(b));
    }
    std::cout << "\nBuffer size: " << buffer.GetSize() << std::endl;
    std::cout << "Print from buffer pointer: " << std::endl;
    auto p = buffer.GetData();
    for (auto i = 0; i < buffer.GetSize(); i++) {
      printf("%02x ", static_cast<unsigned char>(p[i]));
    }

    std::cout << "Print from value.ptr: " << std::endl;
    // value.ptr =
    //     std::shared_ptr<const char[]>(buffer.GetPointer(), [](const char *)
    //     {});
    value.ptr = buffer.GetData();
    for (auto i = 0; i < buffer.GetSize(); i++) {
      printf("%02x ", static_cast<unsigned char>(value.ptr[i]));
    }

    std::cout << "value.ptr == buffer.GetData(): "
              << (value.ptr == buffer.GetData()) << std::endl;
  }

  std::string GetBinaryMol() {
    std::string buffer;
    buffer.resize(value.length - PREFIX_LENGTH);
    if (value.ptr && value.length > PREFIX_LENGTH) {
      memcpy(&buffer[0], &value.ptr[PREFIX_LENGTH],
             value.length - PREFIX_LENGTH);
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
