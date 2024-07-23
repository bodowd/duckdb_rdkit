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
  static constexpr idx_t NUM_ATOMS_BYTES = 4;
  static constexpr idx_t NUM_BONDS_BYTES = 4;
  static constexpr idx_t AMW_BYTES = 4;
  static constexpr idx_t NUM_RINGS_BYTES = 4;
  static constexpr idx_t BMOL_SIZE_BYTES = 4;
  static constexpr idx_t HEADER_SIZE = NUM_ATOMS_BYTES + NUM_BONDS_BYTES +
                                       AMW_BYTES + NUM_RINGS_BYTES +
                                       BMOL_SIZE_BYTES;
  static constexpr idx_t UINT32_MAX_SIZE = NumericLimits<uint32_t>::Maximum();
  uint32_t num_atoms;
  uint32_t num_bonds;
  uint32_t amw;
  uint32_t num_rings;
  uint32_t bmol_size;
  std::string bmol;

  // default constructor for deserialization
  umbra_mol_t() = default;

  umbra_mol_t(uint32_t num_atoms, uint32_t num_bonds, uint32_t amw,
              uint32_t num_rings, const std::string &binary_mol)
      : num_atoms(num_atoms), num_bonds(num_bonds), amw(amw),
        num_rings(num_rings), bmol_size(binary_mol.size()), bmol(binary_mol) {

    if (num_atoms > UINT32_MAX_SIZE || num_bonds > UINT32_MAX_SIZE ||
        amw > UINT32_MAX_SIZE || num_rings > UINT32_MAX_SIZE) {
      throw OutOfRangeException(
          "Cannot support a molecule of this size. There are properties of "
          "this molecule larger than the supported size: '%d'",
          UINT32_MAX_SIZE);
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
