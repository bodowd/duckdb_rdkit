#pragma once
#include "common.hpp"
#include "duckdb/common/exception.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Fingerprints/MACCS.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <bitset>
#include <cstdint>
#include <memory>
#include <string>
#include <sys/types.h>

namespace duckdb_rdkit {

LogicalType Mol();
LogicalType UmbraMol();
void RegisterTypes(DatabaseInstance &instance);

struct umbra_mol_t {

private:
  // fragments and number of times they appear
  // use a vector to keep the order. The order oof the vector will inform the
  // bit order
  // "O" 2 times, is bit 0,
  // "O" 3 times is bit 1, etc...
  std::vector<std::vector<string>> dalke_fp = {{"O", "2", "3", "1", "4", "5"},
                                               {"Ccc", "2", "4"},
                                               {"CCN", "1"},
                                               {"cnc", "1"},
                                               {"cN", "1"},
                                               {"C=O", "1"},
                                               {"CCC", "1"},
                                               {"S", "1"},
                                               {"c1ccccc1", "1", "2"},
                                               {"N", "2", "3", "1"},
                                               {"C=C", "1"},
                                               {"nn", "1"},
                                               {"CO", "2"},
                                               {"Ccn", "1", "2"},
                                               {"CCCCC", "1"},
                                               {"cc(c)c", "1"},
                                               {"CNC", "2"},
                                               {"s", "1"},
                                               {"CC(C)C", "1"},
                                               {"o", "1"},
                                               {"cncnc", "1"},
                                               {"C=N", "1"},
                                               {"CC=O", "2", "3"},
                                               {"Cl", "1"},
                                               {"ccncc", "2"},
                                               {"CCCCCC", "6"},
                                               {"F", "1"},
                                               {"CCOC", "3"},
                                               {"c(cn)n", "1"},
                                               {"C", "9", "6", "1"},
                                               {"CC=C(C)C", "1"},
                                               {"c1ccncc1", "1"},
                                               {"CC(C)N", "1"},
                                               {"CC", "1"},
                                               {"CCC(C)O", "4"},
                                               {"ccc(cc)n", "2"},
                                               {"C1CCCC1", "1"},
                                               {"CNCN", "1"},
                                               {"cncn", "3"},
                                               {"CSC", "1"},
                                               {"CCNCCCN", "1"},
                                               {"CccC", "1"},
                                               {"ccccc(c)c", "3"}};

  void make_dalke_fp() {
    std::unique_ptr<RDKit::ROMol> target(new RDKit::ROMol());
    RDKit::MolPickler::molFromPickle(bmol, *target);
    // std::cout << "Making Dalke FP" << std::endl;
    RDKit::SubstructMatchParameters params;
    params.uniquify = true;
    params.useQueryQueryMatches = false;
    params.recursionPossible = true;
    params.useChirality = false;
    params.maxMatches = 10;
    params.numThreads = 1;

    uint8_t curBit = 0;
    // for each of the dalke fragments, check if it is found
    // in the target molecule, the one that an UmbraMol will be constructed for
    // the dalke fp is the query molecule
    for (const auto &fp : dalke_fp) {
      std::unique_ptr<RDKit::ROMol> dalke_fp_mol;

      try {
        dalke_fp_mol.reset(RDKit::SmilesToMol(fp[0], 0, false));
      } catch (std::exception &e) {
        std::string msg = StringUtil::Format("%s", typeid(e).name());
        throw FatalException(msg);
      }

      auto matchVect = RDKit::SubstructMatch(*target, *dalke_fp_mol, params);

      // if the target has the fp substructure in it at least $number of times
      // it appears, set that bit
      for (auto i = 1; i < fp.size(); i++) {
        if (matchVect.size() >= std::stoi(fp[i])) {
          dalke_bitset.set(curBit);
        }

        curBit++;
      }
    }

    D_ASSERT(curBit == 55);

    ;
  }

public:
  static constexpr idx_t NUM_ATOMS_BYTES = 2;
  static constexpr idx_t NUM_BONDS_BYTES = 2;
  static constexpr idx_t AMW_BYTES = 2;
  static constexpr idx_t NUM_RINGS_BYTES = 2;
  static constexpr idx_t BMOL_SIZE_BYTES = 2;
  static constexpr idx_t DALKE_BIT_VECT_SIZE_BITS = 56;
  static constexpr idx_t DALKE_BIT_VECT_SIZE_BYTES =
      7; // 55 bits in the bitvector, but 8 bytes (56 bits) is the closest in
         // bytes, and easier to do deserialization with bytes alignment
  static constexpr idx_t HEADER_SIZE =
      NUM_ATOMS_BYTES + NUM_BONDS_BYTES + AMW_BYTES + NUM_RINGS_BYTES +
      BMOL_SIZE_BYTES + DALKE_BIT_VECT_SIZE_BYTES;
  static constexpr idx_t UINT16_MAX_SIZE = NumericLimits<uint16_t>::Maximum();
  uint16_t num_atoms;
  uint16_t num_bonds;
  uint16_t amw;
  uint16_t num_rings;
  uint16_t bmol_size;
  std::bitset<DALKE_BIT_VECT_SIZE_BITS> dalke_bitset;
  std::string bmol;

  // default constructor for deserialization
  umbra_mol_t() = default;
  // umbra_mol_t(uint16_t num_atoms, uint16_t num_bonds, uint16_t amw,
  //             uint16_t num_rings, const std::string &binary_mol,
  //             const RDKit::ROMol &mol, bool dalke_fp)
  //     : num_atoms(num_atoms), num_bonds(num_bonds), amw(amw),
  //       num_rings(num_rings), bmol_size(binary_mol.size()), bmol(binary_mol)
  //       {
  //
  //   if (num_atoms > UINT16_MAX_SIZE || num_bonds > UINT16_MAX_SIZE ||
  //       amw > UINT16_MAX_SIZE || num_rings > UINT16_MAX_SIZE) {
  //     throw OutOfRangeException(
  //         "Cannot support a molecule of this size. There are properties of "
  //         "this molecule larger than the supported size: '%d'",
  //         UINT16_MAX_SIZE);
  //   }
  //
  //   if (dalke_fp) {
  //     GenerateDalkeFP();
  //   }
  // }

  umbra_mol_t(uint16_t num_atoms, uint16_t num_bonds, uint16_t amw,
              uint16_t num_rings, const std::string &binary_mol,
              const RDKit::ROMol &mol)
      : num_atoms(num_atoms), num_bonds(num_bonds), amw(amw),
        num_rings(num_rings), bmol_size(binary_mol.size()), bmol(binary_mol) {

    if (num_atoms > UINT16_MAX_SIZE || num_bonds > UINT16_MAX_SIZE ||
        amw > UINT16_MAX_SIZE || num_rings > UINT16_MAX_SIZE) {
      throw OutOfRangeException(
          "Cannot support a molecule of this size. There are properties of "
          "this molecule larger than the supported size: '%d'",
          UINT16_MAX_SIZE);
    }
  }

  void GenerateDalkeFP() { make_dalke_fp(); }

  friend std::ostream &operator<<(std::ostream &out,
                                  const umbra_mol_t &umbra_mol) {
    out << "num_atoms: " << umbra_mol.num_atoms << '\n';
    out << "num_bonds: " << umbra_mol.num_bonds << '\n';
    out << "amw: " << umbra_mol.amw << '\n';
    out << "num_rings: " << umbra_mol.num_rings << '\n';
    out << "bmol_size: " << umbra_mol.bmol_size << '\n';
    out << "dalke_bit_vect" << '\n';
    out << umbra_mol.dalke_bitset.to_string() << '\n';
    out << "\nbmol: " << '\n';
    for (char byte : umbra_mol.bmol) {
      printf("%02x ", static_cast<unsigned char>(byte));
    }
    return out;
  }
};

} // namespace duckdb_rdkit
