#pragma once
#include "common.hpp"
#include "duckdb/common/exception.hpp"
#include <cstddef>
#include <cstdint>
#include <sys/types.h>

namespace duckdb_rdkit {

LogicalType Mol();
LogicalType UmbraMol();
void RegisterTypes(DatabaseInstance &instance);

struct umbra_mol_t {

public:
  uint32_t prefix;
  std::string bmol;

  umbra_mol_t(uint32_t num_atoms, uint32_t num_bonds, uint32_t amw,
              uint32_t num_rings, const std::string &binary_mol)
      : bmol(binary_mol) {

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
    // shift 25 bits to the left to pack up to 32 bit (will fill bits 25 to 32)
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
    // 0x7FF is 2047 is 0111_1110_1000	which will set the lowest 11 bits
    // 32 - 7 - 6 - 3 - 11 = 5
    // Shift to the left so that the last bit is at the 5th place
    prefix |= (amw & 0x7FF) << 5;
  }

  friend std::ostream &operator<<(std::ostream &out,
                                  const umbra_mol_t &umbra_mol) {
    // 0x07F is 0111 1111
    // the number & 0x07F will get all lower 7 bits that are set
    auto num_atoms = (umbra_mol.prefix >> 25) & 0x7F;
    auto num_bonds = (umbra_mol.prefix >> 19) & 0x3F;
    auto num_rings = (umbra_mol.prefix >> 16) & 0x07;
    auto amw = (umbra_mol.prefix >> 5) & 0x7FF;
    out << "num_atoms: " << num_atoms << '\n';
    out << "num_bonds: " << num_bonds << '\n';
    out << "num_rings: " << num_rings << '\n';
    out << "amw: " << amw << '\n';
    out << "bmol: " << '\n';
    for (char byte : umbra_mol.bmol) {
      printf("%02x ", static_cast<unsigned char>(byte));
    }
    return out;
  }
};

} // namespace duckdb_rdkit
