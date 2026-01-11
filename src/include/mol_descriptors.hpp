#pragma once
#include "common.hpp"
#include "duckdb/main/connection.hpp"

namespace duckdb_rdkit {
void RegisterDescriptorFunctions(DatabaseInstance &instance);
void mol_registration_hash(DataChunk &args, ExpressionState &state, Vector &result);
}
