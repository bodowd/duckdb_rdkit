#pragma once

#include "duckdb/main/extension/extension_loader.hpp"
#include <GraphMol/GraphMol.h>
#include <memory>

namespace duckdb {

// these functions are used in other parts of the extension, for example in
// casts
std::unique_ptr<RDKit::ROMol> rdkit_mol_from_smiles(std::string s);
std::string rdkit_mol_to_binary_mol(const RDKit::ROMol mol);
std::unique_ptr<RDKit::ROMol> rdkit_binary_mol_to_mol(std::string bmol);
std::string rdkit_mol_to_smiles(RDKit::ROMol mol);

void RegisterFormatFunctions(ExtensionLoader &loader);
} // namespace duckdb
