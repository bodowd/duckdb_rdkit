#pragma once
#include "common.hpp"
#include "duckdb/function/function.hpp"
#include <GraphMol/GraphMol.h>

namespace duckdb_rdkit {

// these functions are used in other parts of the extension, for example in
// casts
std::unique_ptr<RDKit::ROMol> rdkit_mol_from_smiles(std::string smiles);
std::string rdkit_mol_to_binary_mol(const RDKit::ROMol mol);
std::unique_ptr<RDKit::ROMol> rdkit_binary_mol_to_mol(std::string bmol);
std::string rdkit_mol_to_smiles(RDKit::ROMol mol);
bool mol_cmp(const RDKit::ROMol &m1, const RDKit::ROMol &m2);

struct BindData : public FunctionData {
public:
  BindData(const string &smiles);
  // BindData(const string &smiles, unique_ptr<RDKit::ROMol> mol);
  BindData(const BindData &other);

  bool Equals(const FunctionData &other_p) const override;
  unique_ptr<FunctionData> Copy() const override;

  string smiles;
  unique_ptr<RDKit::ROMol> mol;

  void InitMol();
};

//! Binds a default calendar object for use by the function
static unique_ptr<FunctionData>
Bind(ClientContext &context, ScalarFunction &bound_function,
     vector<duckdb::unique_ptr<Expression>> &arguments);

void RegisterFormatFunctions(DatabaseInstance &instance);
} // namespace duckdb_rdkit
