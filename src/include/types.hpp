#pragma once
#include "common.hpp"
#include "duckdb/common/exception.hpp"
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

  umbra_mol_t(uint num_atoms, uint num_bonds, uint amw, uint num_rings,
              const std::string &binary_mol)
      : bmol(binary_mol) {
    std::cout << "CONSTRUCTOR FOR umbra_mol_t" << std::endl;
    std::cout << "NUM_ATOMS IN INPUT: " << num_atoms << std::endl;
    // cap the count if it is larger than the number of bits it supports
    // number of bits for each count supports the 99 percentile of
    // values in chembl
    if (num_atoms > 128) {
      num_atoms = 128;
    }

    if (num_bonds > 64) {
      num_bonds = 64;
    }

    if (amw > 2048) {
      amw = 2048;
    }

    if (num_rings > 8) {
      num_rings = 8;
    }

    // 0x7F is 0111 1111 which sets the first 7 bits of a number
    // (num_atoms & 0x7F) creates a mask, keeping only the first 7 bits
    // and higher order bits are zeroed out
    // Shift 20 to the left in order to shift to the 20th bit
    // The end of the number will be at 27th bit
    prefix |= (num_atoms & 0x7F) << 20;
    // 0x3F is 0011 1111 to set the first 6 bits
    // apply mask to keep only the first 6 bits
    // shift 14 bits to the left to put it at the 14th bit
    prefix |= (num_bonds & 0x3F) << 14;
    prefix |= (num_rings & 0x07) << 11;
    // 11 bits for amw, doesn't need to shift
    prefix |= (amw & 0x7FF);

    std::bitset<32> p(prefix);
    std::cout << "prefix: " << p << '\n' << std::endl;
  }

  friend std::ostream &operator<<(std::ostream &out,
                                  const umbra_mol_t &umbra_mol) {
    // 0x07F is 0111 1111
    // the number & 0x07F will get all lower 7 bits that are set
    auto num_atoms = (umbra_mol.prefix >> 20) & 0x7F;
    auto num_bonds = (umbra_mol.prefix >> 14) & 0x3F;
    auto num_rings = (umbra_mol.prefix >> 11) & 0x07;
    auto amw = umbra_mol.prefix & 0x7FF;
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
