#pragma once
#include "common.hpp"
#include "duckdb/common/exception.hpp"
#include <cstddef>
#include <cstdint>
#include <sys/types.h>

namespace duckdb_rdkit {

LogicalType Mol();
void RegisterTypes(DatabaseInstance &instance);
} // namespace duckdb_rdkit
