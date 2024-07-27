#pragma once
#include "common.hpp"
#include "duckdb/common/exception.hpp"
#include <GraphMol/Fingerprints/MACCS.h>
#include <cstdint>
#include <memory>
#include <sys/types.h>

namespace duckdb_rdkit {

LogicalType Mol();
LogicalType UmbraMol();
void RegisterTypes(DatabaseInstance &instance);

struct umbra_mol_t {

public:
  static constexpr idx_t NUM_ATOMS_BYTES = 2;
  static constexpr idx_t NUM_BONDS_BYTES = 2;
  static constexpr idx_t AMW_BYTES = 2;
  static constexpr idx_t NUM_RINGS_BYTES = 2;
  static constexpr idx_t BMOL_SIZE_BYTES = 2;
  static constexpr idx_t MACCS_SIZE_BYTES = 2;
  static constexpr idx_t HEADER_SIZE = NUM_ATOMS_BYTES + NUM_BONDS_BYTES +
                                       AMW_BYTES + NUM_RINGS_BYTES +
                                       BMOL_SIZE_BYTES + MACCS_SIZE_BYTES;
  static constexpr idx_t UINT16_MAX_SIZE = NumericLimits<uint16_t>::Maximum();
  uint16_t num_atoms;
  uint16_t num_bonds;
  uint16_t amw;
  uint16_t num_rings;
  uint16_t bmol_size;
  char maccs[MACCS_SIZE_BYTES]; // 167 bits in RDKit is a little less than 21
                                // bytes
  std::string bmol;

  // default constructor for deserialization
  umbra_mol_t() = default;

  umbra_mol_t(uint16_t num_atoms, uint16_t num_bonds, uint16_t amw,
              uint16_t num_rings,
              std::unique_ptr<ExplicitBitVect> maccs_bit_vect,
              const std::string &binary_mol)
      : num_atoms(num_atoms), num_bonds(num_bonds), amw(amw),
        num_rings(num_rings), bmol_size(binary_mol.size()), bmol(binary_mol) {

    // zero out the maccs char array
    std::memset(maccs, 0, MACCS_SIZE_BYTES);
    std::memcpy(maccs, &maccs_bit_vect, MACCS_SIZE_BYTES);
    std::cout << "maccs in constructor: " << std::endl;
    for (auto i = 0; i < MACCS_SIZE_BYTES; i++) {
      printf("%02x ", maccs[i]);
    }

    if (num_atoms > UINT16_MAX_SIZE || num_bonds > UINT16_MAX_SIZE ||
        amw > UINT16_MAX_SIZE || num_rings > UINT16_MAX_SIZE) {
      throw OutOfRangeException(
          "Cannot support a molecule of this size. There are properties of "
          "this molecule larger than the supported size: '%d'",
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
