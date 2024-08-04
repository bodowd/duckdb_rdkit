#pragma once
#include "common.hpp"
#include "duckdb/common/typedefs.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <cstdint>
#include <cstring>
#include <memory>

namespace duckdb {

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

  umbra_mol_t(uint32_t num_atoms, uint32_t num_bonds, uint32_t amw,
              uint32_t num_rings, const std::string &binary_mol) {
    // the length of umbra_mol_t is the prefix and binary mol
    // because all of that will be stored in the ptr
    // This makes the implementation more similar to the string_t, which
    // copies only the first $PREFIX_LENGTH bytes from the string into
    // the prefix, but still contains the whole thing in the ptr
    // See Finalize function in string_type.hpp
    //
    // Since the PREFIX is calculated in the constructor from the binary mol
    // and is not provided by the caller. Therefore, we need to add
    // PREFIX_LENGTH to binary_mol
    //
    // GetSize() will use the length field in the inlined struct
    value.inlined.length = PREFIX_LENGTH + binary_mol.size();
    D_ASSERT(data || GetSize() == 0);
    if (IsInlined()) {
    } else {
      /* The prefix has the following format, from highest order bits
       * to the lowest order bits
       * num atoms 7 bits
       * num bonds 6 bits
       * num rings 3 bits
       * amw 11 bits
       * total 27 bits
       * pack into 32 bit unsigned int
       * 5 bits left over for something else if needed
       */

      // zero initialize the prefix
      // otherwise, prefix will contain whatever data was previously in
      // this memory location!
      memset(value.pointer.prefix, 0, PREFIX_LENGTH);

      uint32_t prefix;
      prefix = 0;
      // std::cout << "CONSTRUCTOR FOR umbra_mol_t" << std::endl;
      // std::cout << "NUM_ATOMS IN INPUT: " << num_atoms << std::endl;
      // std::cout << "NUM_BONDS IN INPUT: " << num_bonds << std::endl;
      // std::cout << "NUM_RINGS IN INPUT: " << num_rings << std::endl;
      // std::cout << "AMW IN INPUT: " << amw << std::endl;
      // cap the count if it is larger than the number of bits it supports
      // number of bits for each count supports the 99 percentile of
      // values in chembl
      if (num_atoms >= 127) {
        num_atoms = 127;
      }

      if (num_bonds >= 63) {
        num_bonds = 63;
      }

      if (num_rings >= 7) {
        num_rings = 7;
      }

      if (amw >= 2047) {
        amw = 2047;
      }

      // 0x7F is 127 is 0111 1111 which sets the first 7 bits of a number
      // (num_atoms & 0x7F) creates a mask, keeping only the first 7 bits
      // and higher order bits are zeroed out
      // shift 25 bits to the left to pack up to 32 bit (will fill bits 25 to
      // 32)
      prefix |= (num_atoms & 0x7F) << 25;
      // 0x3F is 63 is 0011 1111 to set the first 6 bits
      // apply mask to keep only the first 6 bits
      // shift 19 bits to the left
      // because first 7 bits are now taken
      // This value should take the next 6 bits.
      // that means 13 bits will be occupied
      // So then set the lowest bit of the num_bonds value to the 19th
      // bit because
      // 32 bits - 13 bits = 19 bits
      prefix |= (num_bonds & 0x3F) << 19;
      // 0x07 is 0111 is 7 which will set the lower 3 bits
      // 32 bits - 7 bits for num_atoms - 6 bits for num_bonds - 3 bits for
      // num_rings = 16 So this should be shifted so that the last bit of
      // num_rings is at the 16th place
      prefix |= (num_rings & 0x07) << 16;
      // 0x7FF is 2047 is 0111_1110_1000	which will set the lowest 11
      // bits 32 - 7 - 6 - 3 - 11 = 5 Shift to the left so that the last bit is
      // at the 5th place
      prefix |= (amw & 0x7FF) << 5;

      // TODO:  dalke_fp

      uint64_t dalke_fp;
      dalke_fp = 0;

      //
      // Only the highest 27 bits are part of the prefix
      // At this point, the prefix looks like this:
      // |-------------- 32 bits ---------------|
      // |--27 bits count prefix--|
      //
      // A uint64_t will be used to construct the dalke_fp, which is has a
      // length of 55 bits.
      // |----------------------64 bits ----------------------|
      // |------- 55 bits dalke_fp -------------|
      //
      // Then, the first 5 bits will be copied to the lower 5 bits of the 32
      // already allocated above.
      // |-------------- 32 bits --------------------|
      // |--27 bits count prefix--|--5 bit dalke_fp--|
      // and
      // The remaining highest order 50 bits will be copied to the prefix as
      // well
      //
      // |-------- 88 bits prefix -------|
      // |- 82 bits counts + dalke_fp -|
      // |----27-----|--5--|----50 ----|
      //
      //
      // In total, the prefix will then be 82 bits (27 bit count prefix + 55 bit
      // dalke_fp), which goes into 11 bytes (88 bits). The lowest order 6 bits
      // of the prefix do not represent anything.

      // The value.pointer.ptr field should contain the prefix+binary mol
      // Therefore, we copy the prefix and binary mol into a std::string buffer
      // in order to combine the data, and then we copy just the first
      // PREFIX_LENGTH bytes of the combined data into the prefix, and
      // then copy the entire combined data into the ptr field
      // Rather than copying prefix just into the value.pointer.prefix field
      // and then the ptr only containing the binary molecule,
      // doing it this way is more consistent with how the string_t
      // implementation is done, and we want to make sure the implemented
      // functions behave like the string_t when the rest of duckdb thinks
      // umbra_mol_t is a string_t
      // First, combine the prefix and binary mol
      // FIXME: Only 27 bits of the prefix is calculated, and the remaining
      // 55 bits will be random numbers that change with each run because
      // dalke_fp is not implemented yet

      // std::string buffer;
      // buffer.resize(PREFIX_BYTES + binary_mol.size());
      // idx_t offset = 0;
      // // copy the first 32 bits of the prefix -- the counts
      // memcpy(&buffer[offset], &prefix, COUNT_PREFIX_BYTES);
      // offset += COUNT_PREFIX_BYTES;
      // // copy the dalke_fp
      // memcpy(&buffer[offset], &dalke_fp, DALKE_FP_PREFIX_BYTES);
      // offset += DALKE_FP_PREFIX_BYTES;
      // // copy the binary mol
      // memcpy(&buffer[offset], &binary_mol[0], binary_mol.size());
      // std::cout << "after memcpy first" << std::endl;
      // // now copy just the first PREFIX_LENGTH bytes to prefix
      // // and the whole combined data to ptr
      // memcpy(value.pointer.prefix, &buffer.data()[0], PREFIX_LENGTH);
      // std::cout << "after memcpy second" << std::endl;
      // value.pointer.ptr = (char *)buffer.data();
      // std::cout << "done" << std::endl;

      std::string buffer;
      buffer.reserve(PREFIX_BYTES + binary_mol.size());
      buffer.append(reinterpret_cast<const char *>(&prefix),
                    COUNT_PREFIX_BYTES);
      buffer.append(reinterpret_cast<const char *>(&dalke_fp),
                    DALKE_FP_PREFIX_BYTES);
      buffer.append(binary_mol);
      // now copy just the first PREFIX_LENGTH bytes to prefix
      // and the whole combined data to ptr
      memcpy(value.pointer.prefix, &buffer.data()[0], PREFIX_LENGTH);
      std::cout << "after memcpy second" << std::endl;
      for (char b : buffer) {
        printf("%02x ", static_cast<unsigned char>(b));
      }
      value.pointer.ptr = buffer.data();

      std::cout << "after setting pointer" << std::endl;

      for (idx_t i = 0; i < PREFIX_LENGTH + binary_mol.size(); i++) {
        printf("%02x ", static_cast<unsigned char>(value.pointer.ptr[i]));
      }

      std::cout << "done" << std::endl;
    }
  }

  // this is used for converting a string_t which is expected to be an
  // umbra_mol_t when it is
  // received by a function processing string vectors
  // for example, with `umbra_is_exact_match``
  umbra_mol_t(const char *data, uint32_t len) {
    value.inlined.length = len;
    // zero out the prefix and copy the first PREFIX_LENGTH bytes to
    // the prefix field
    memset(value.pointer.prefix, 0, PREFIX_LENGTH);
    memcpy(value.pointer.prefix, data, PREFIX_LENGTH);

    // copy from offeset = PREFIX_LENGTH in the data to the end,
    // which is calculated by GetBmolSize, to the pointer
    memcpy(value.pointer.ptr, &data[PREFIX_LENGTH], GetBmolSize());
  }

  umbra_mol_t(const string &
                  value) // NOLINT: Allow implicit conversion from `const char*`
      : umbra_mol_t(value.c_str(), UnsafeNumericCast<uint32_t>(value.size())) {}

  bool IsInlined() const { return GetSize() <= INLINE_LENGTH; }

  // returns just the binary molecule, not the prefix
  const char *GetBinaryMol() const { return GetData(); }

  const char *GetData() const {
    std::cout << "Prefix" << std::endl;
    for (auto i = 0; i < PREFIX_LENGTH; i++) {
      printf("%02x ", value.pointer.prefix[i]);
    }

    std::cout << "Pointer: " << std::endl;
    for (auto i = 0; i < GetSize(); i++) {
      printf("%02x ", static_cast<unsigned char>(value.pointer.ptr[i]));
    }

    return IsInlined() ? const_char_ptr_cast(value.inlined.inlined)
                       : value.pointer.ptr;
  }
  const char *GetDataUnsafe() const { return GetData(); }

  char *GetDataWriteable() const {
    return IsInlined() ? (char *)value.inlined.inlined
                       : value.pointer.ptr; // NOLINT
  }

  const char *GetPrefix() const {
    // I don't understand, but somehow value.inlined.inlined also gets populated
    // even if it is a big int
    return value.inlined.inlined;
  }

  char *GetPrefixWriteable() { return value.inlined.inlined; }

  // This is PREFIX_LENGTH + binary mol size
  idx_t GetSize() const { return value.inlined.length; }

  idx_t GetBmolSize() const { return GetSize() - PREFIX_LENGTH; }

  bool Empty() const { return value.inlined.length == 0; }

  string GetString() const { return string(GetData(), GetSize()); }

  explicit operator string() const { return GetString(); }

  char *GetPointer() const {
    D_ASSERT(!IsInlined());
    return value.pointer.ptr;
  }

  void SetPointer(char *new_ptr) {
    D_ASSERT(!IsInlined());
    value.pointer.ptr = new_ptr;
  }

  void Finalize() {
    auto dataptr = GetData();
    memcpy(value.pointer.prefix, dataptr, PREFIX_LENGTH);
  }

  void Verify() const;
  void VerifyUTF8() const;
  void VerifyCharacters() const;
  void VerifyNull() const;

private:
  union {
    struct {
      // this is the binary molecule size
      // the prefix will be calculated in the umbra_mol_t struct constructor
      uint32_t length;
      char prefix[PREFIX_BYTES];
      char *ptr;
    } pointer;
    struct {
      uint32_t length;
      char inlined[12];
    } inlined;
  } value;
};

} // namespace duckdb
