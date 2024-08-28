#pragma once
#include "umbra_mol.hpp"
#include "common.hpp"
#include "duckdb/common/exception.hpp"
#include "mol_formats.hpp"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <cstdint>
#include <string>

namespace duckdb_rdkit {

// These are fragments and the number of times they appear in a molecule.
// Use a vector to keep the order. The order of the vector will inform the
// bit order.
// "O" 2 times, is bit 0,
// "O" 3 times is bit 1, etc...
// This was calculated by Andrew Dalke:
// http://www.dalkescientific.com/writings/diary/archive/2012/06/11/optimizing_substructure_keys.html
// And Greg Landrum tested out the 55 bits I use here.
// https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg02078.html
//
// This is used for a substructure filter and placed in the prefix of a
// Umbra-mol so that short-circuiting can take place and save on computation
// cost when deserializing and a full substructure search is unnecessary.
std::vector<std::vector<string>> dalke_counts = {{"O", "2", "3", "1", "4", "5"},
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

uint64_t make_dalke_fp(const RDKit::ROMol &mol) {
  std::bitset<64> bs;
  RDKit::SubstructMatchParameters params;
  params.uniquify = true;
  params.useQueryQueryMatches = false;
  params.recursionPossible = true;
  params.useChirality = false;
  params.maxMatches = 10;
  params.numThreads = 1;

  uint8_t curBit = 0;
  // for each of the dalke fragments, check if it is found
  // in the target molecule (the one that an UmbraMol will be constructed for)
  // the dalke fragment is the "query" molecule in the SubstructMatch
  // function
  for (const auto &fp : dalke_counts) {
    std::unique_ptr<RDKit::ROMol> dalke_fp_mol;

    try {
      dalke_fp_mol.reset(RDKit::SmilesToMol(fp[0], 0, false));
    } catch (std::exception &e) {
      std::string msg = StringUtil::Format("%s", typeid(e).name());
      throw InvalidInputException(msg);
    }

    auto matchVect = RDKit::SubstructMatch(mol, *dalke_fp_mol, params);

    // if the target has the fp substructure in it at least $NUMBER of times
    // it appears, set that bit
    for (auto i = 1; i < fp.size(); i++) {
      if (matchVect.size() >= std::stoi(fp[i])) {
        bs.set(curBit);
      }

      curBit++;
    }
  }
  auto k = bs.to_ullong();
  D_ASSERT(curBit == 55);
  return k;
}

// "Umbra-mol" has more than just the binary molecule
// There is a prefix in front of the binary molecule, inspired by
// Umbra-style strings
std::string get_umbra_mol_string(const RDKit::ROMol &mol) {
  auto binary_mol = rdkit_mol_to_binary_mol(mol);
  size_t total_size = umbra_mol_t::PREFIX_BYTES + binary_mol.size();

  // Add the prefix in front of the binary molecule object
  //
  // The prefix consists of two parts: some simple counts used for filtering
  // records in exact match and a bit vector that marks the presence of
  // certain substructures for filtering records in substructure match
  // searches

  // The counts prefix consists of
  // number of atoms, 7 bits
  // number of bonds, 6 bits
  // number of rings, 3 bits
  // amw, 11 bits
  // In total that is 27 bits, which is placed into a uint32_t.
  // There are 5 bits left over, but this makes alignment a little easier
  // particularly when the dalke fingerprint vector is included
  uint32_t num_atoms = mol.getNumAtoms();
  uint32_t num_bonds = mol.getNumBonds();
  uint32_t amw = RDKit::Descriptors::calcAMW(mol);
  uint32_t num_rings = mol.getRingInfo()->numRings();

  // initialize with zero otherwise prefix will have data from
  // previous calls of this function
  uint32_t prefix = 0;

  // Cap the count if it is larger than the number of bits it supports.
  // The number of bits to use was determined by an analysis of these properties
  // of molecules in chembl.
  // The number of bits for each count supports the 99th percentile of
  // values in chembl.
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

  uint64_t dalke_fp = make_dalke_fp(mol);

  // remember to keep endianess in mind if you print things out.
  // little endian on my machine
  std::string buffer;
  buffer.reserve(total_size);
  buffer.append(reinterpret_cast<const char *>(&prefix),
                umbra_mol_t::COUNT_PREFIX_BYTES);
  buffer.append(reinterpret_cast<const char *>(&dalke_fp),
                umbra_mol_t::DALKE_FP_PREFIX_BYTES);
  buffer.append(binary_mol);

  return buffer;
}
} // namespace duckdb_rdkit
