#pragma once
#include "common.hpp"

namespace duckdb_rdkit {

void VarcharToMol(Vector &source, Vector &result, idx_t count);
bool VarcharToMolCast(Vector &source, Vector &result, idx_t count,
                      CastParameters &parameters);
void MolToVarchar(Vector &source, Vector &result, idx_t count);
bool MolToVarcharCast(Vector &source, Vector &result, idx_t count,
                      CastParameters &parameters);
void RegisterCasts(DatabaseInstance &instance);

} // namespace duckdb_rdkit
