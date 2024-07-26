#pragma once
#include "common.hpp"
#include "duckdb/common/exception.hpp"
#include <cstddef>
#include <sys/types.h>

namespace duckdb_rdkit {

LogicalType Mol();
LogicalType UmbraMol();
void RegisterTypes(DatabaseInstance &instance);

struct umbra_mol_t {

public:
  char prefix[8];
  std::string bmol;

  static constexpr idx_t PREFIX_SIZE = 8;

  // default constructor for deserialization
  umbra_mol_t() = default;

  umbra_mol_t(const char *canonical_smiles, const std::string &binary_mol)
      : bmol(binary_mol) {

    std::size_t len = std::strlen(canonical_smiles);
    // zero out the prefix
    std::memset(prefix, 0, PREFIX_SIZE);

    // for (auto i = 0; i < len; i++) {
    //   printf("%02x ", static_cast<unsigned char>(canonical_smiles[i]));
    // }
    std::memcpy(prefix, canonical_smiles, PREFIX_SIZE);
    if (len < PREFIX_SIZE) {
      // fill the prefix with zeros if the len of the smiles string is <
      // PREFIX_SIZE
      std::memset(prefix + len, 0, PREFIX_SIZE - len);
    }
  }

  friend std::ostream &operator<<(std::ostream &out,
                                  const umbra_mol_t &umbra_mol) {
    out << "prefix: " << '\n';
    for (char i : umbra_mol.prefix) {
      printf("%02x ", static_cast<unsigned char>(i));
    }
    out << "\nbmol: " << '\n';
    for (char byte : umbra_mol.bmol) {
      printf("%02x ", static_cast<unsigned char>(byte));
    }
    return out;
  }
};

} // namespace duckdb_rdkit
