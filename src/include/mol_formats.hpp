#pragma once
#include "common.hpp"
#include "types.hpp"
#include <GraphMol/GraphMol.h>

namespace duckdb_rdkit {

// these functions are used in other parts of the extension, for example in
// casts
std::unique_ptr<RDKit::ROMol> rdkit_mol_from_smiles(std::string s);
std::string rdkit_mol_to_binary_mol(const RDKit::ROMol mol);
std::unique_ptr<RDKit::ROMol> rdkit_binary_mol_to_mol(std::string bmol);
std::string rdkit_mol_to_smiles(RDKit::ROMol mol);
std::string serialize_umbra_mol(umbra_mol_t umbra_mol);
uint32_t extract_prefix_from_umbra_mol(string_t buffer);
std::string extract_bmol_from_umbra_mol(string_t buffer);

void RegisterFormatFunctions(DatabaseInstance &instance);
} // namespace duckdb_rdkit
