#pragma once
#include "common.hpp"
#include "duckdb/common/assert.hpp"
#include "duckdb/common/constants.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/serializer/binary_serializer.hpp"
#include "duckdb/common/serializer/memory_stream.hpp"
#include "duckdb/common/typedefs.hpp"
#include <cstdint>
#include <memory>
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
  size_t bmol_size;
  std::string bmol;

  // default constructor for deserialization
  umbra_mol_t() = default;

  umbra_mol_t(uint32_t num_atoms, uint32_t num_bonds, uint32_t amw,
              uint32_t num_rings, const std::string &binary_mol)
      : num_atoms(num_atoms), num_bonds(num_bonds), amw(amw),
        num_rings(num_rings), bmol_size(binary_mol.size()), bmol(binary_mol) {
    // std::cout << "making mol_t" << std::endl;
    // std::cout << "num_atoms: " << num_atoms << std::endl;
    // std::cout << "num_bonds: " << num_bonds << std::endl;
    // std::cout << "amw: " << amw << std::endl;
    // std::cout << "num_rings: " << num_rings << std::endl;
    // std::cout << "binary_mol: " << std::endl;
    // for (char byte : binary_mol) {
    //   printf("%02x ", static_cast<unsigned char>(byte));
    // }
    // for (size_t i = 0; i < bmol_size; ++i) {
    //   auto a = bmol;
    //   printf("%02x ", static_cast<unsigned char>(a[i]));
    // }

    if (num_atoms > UINT32_MAX_SIZE || num_bonds > UINT32_MAX_SIZE ||
        amw > UINT32_MAX_SIZE || num_rings > UINT32_MAX_SIZE) {
      throw OutOfRangeException(
          "Cannot support a molecule of this size. There are properties of "
          "this molecule larger than the supported size: '%d'",
          UINT32_MAX_SIZE);
    }
  }
};

} // namespace duckdb_rdkit
