#pragma once

#include "duckdb/common/types/vector.hpp"
#include "duckdb/function/cast/default_casts.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

void VarcharToMol(Vector &source, Vector &result, idx_t count);
bool VarcharToMolCast(Vector &source, Vector &result, idx_t count,
                      CastParameters &parameters);
void MolToVarchar(Vector &source, Vector &result, idx_t count);
bool MolToVarcharCast(Vector &source, Vector &result, idx_t count,
                      CastParameters &parameters);
void RegisterCasts(ExtensionLoader &loader);

} // namespace duckdb
