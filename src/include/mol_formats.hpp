#pragma once
#include "common.hpp"
#include <GraphMol/GraphMol.h>

namespace duckdb_rdkit {

// these functions are used in other parts of the extension, for example in
// casts
RDKit::ROMol *rdkit_mol_from_smiles(std::string s);
std::string rdkit_mol_to_binary_mol(const RDKit::ROMol mol);
RDKit::ROMol rdkit_binary_mol_to_mol(std::string bmol);
std::string rdkit_mol_to_smiles(RDKit::ROMol mol);

void RegisterFunctions(DatabaseInstance &instance);
} // namespace duckdb_rdkit
