#pragma once
#include "common.hpp"

namespace duckdb_rdkit {

LogicalType Mol();
void RegisterTypes(DatabaseInstance &instance);

} // namespace duckdb_rdkit
