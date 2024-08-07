#pragma once
#include "umbra_mol.hpp"
#include "common.hpp"
#include <cstdint>
#include <string>

namespace duckdb_rdkit {
std::string get_umbra_mol_string(uint32_t num_atoms, uint32_t num_bonds,
                                 uint32_t amw, uint32_t num_rings,
                                 const std::string &binary_mol) {
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

  std::string buffer;
  buffer.reserve(umbra_mol_t::PREFIX_BYTES + binary_mol.size());
  buffer.append(reinterpret_cast<const char *>(&prefix),
                umbra_mol_t::COUNT_PREFIX_BYTES);
  buffer.append(reinterpret_cast<const char *>(&dalke_fp),
                umbra_mol_t::DALKE_FP_PREFIX_BYTES);
  buffer.append(binary_mol);

  // for (char b : buffer) {
  //   printf("%02x ", static_cast<unsigned char>(b));
  // }

  return buffer;
}
} // namespace duckdb_rdkit
