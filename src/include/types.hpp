#pragma once
#include "common.hpp"

namespace duckdb_rdkit {

LogicalType Mol();

void Register(DatabaseInstance &db);

} // namespace duckdb_rdkit
