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
  static constexpr idx_t NUM_ATOMS_BYTES = 1;
  static constexpr idx_t NUM_BONDS_BYTES = 1;
  static constexpr idx_t AMW_BYTES = 2;
  static constexpr idx_t NUM_RINGS_BYTES = 1;
  static constexpr idx_t BMOL_SIZE_BYTES = 2;
  static constexpr idx_t HEADER_SIZE = NUM_ATOMS_BYTES + NUM_BONDS_BYTES +
                                       AMW_BYTES + NUM_RINGS_BYTES +
                                       BMOL_SIZE_BYTES;
  static constexpr idx_t UINT16_MAX_SIZE = NumericLimits<uint16_t>::Maximum();
  static constexpr idx_t UINT8_MAX_SIZE = NumericLimits<uint8_t>::Maximum();

  // Experimented with 1 byte values 7 byte header
  // 2 bytes also works fine to give 10 byte header
  uint8_t num_atoms;
  uint8_t num_bonds;
  uint16_t amw;
  uint8_t num_rings;
  uint16_t bmol_size;
  std::string bmol;

  // default constructor for deserialization
  umbra_mol_t() = default;

  umbra_mol_t(uint8_t num_atoms, uint8_t num_bonds, uint16_t amw,
              uint8_t num_rings, const std::string &binary_mol)
      : num_atoms(num_atoms), num_bonds(num_bonds), amw(amw),
        num_rings(num_rings), bmol_size(binary_mol.size()), bmol(binary_mol) {

    if (num_atoms > UINT8_MAX_SIZE) {
      num_atoms = UINT8_MAX_SIZE;
    }

    if (num_bonds > UINT8_MAX_SIZE) {
      num_bonds = UINT8_MAX_SIZE;
    }

    if (amw > UINT16_MAX_SIZE) {
      amw = UINT16_MAX_SIZE;
    }

    if (num_rings > UINT8_MAX_SIZE) {
      num_rings = UINT8_MAX_SIZE;
    }

    if (bmol_size > UINT16_MAX_SIZE) {
      throw OutOfRangeException(
          "Cannot support a molecule of this size."
          "this molecule object is larger than the max supported size: '%d'",
          UINT16_MAX_SIZE);
    }
  }

  friend std::ostream &operator<<(std::ostream &out,
                                  const umbra_mol_t &umbra_mol) {
    out << "num_atoms: " << umbra_mol.num_atoms << '\n';
    out << "num_bonds: " << umbra_mol.num_bonds << '\n';
    out << "amw: " << umbra_mol.amw << '\n';
    out << "num_rings: " << umbra_mol.num_rings << '\n';
    out << "bmol_size: " << umbra_mol.bmol_size << '\n';
    out << "bmol: " << '\n';
    for (char byte : umbra_mol.bmol) {
      printf("%02x ", static_cast<unsigned char>(byte));
    }
    return out;
  }
};

} // namespace duckdb_rdkit
