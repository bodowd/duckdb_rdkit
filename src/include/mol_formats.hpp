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
umbra_mol_t deserialize_umbra_mol(std::string buffer);
std::string deserialize_umbra_mol_header(std::string buffer);
std::string deserialize_umbra_mol_bmol(std::string buffer);
std::string serialize_umbra_mol(umbra_mol_t umbra_mol);

void RegisterFormatFunctions(DatabaseInstance &instance);
} // namespace duckdb_rdkit
