#pragma once

#include "duckdb.hpp"
#include "duckdb/common/helper.hpp"
#include "duckdb/main/extension_util.hpp"

// including common.hpp into the other files makes it so that
// it is not necessary to put duckdb::FUNCTION. Brings in the namespace
using namespace duckdb;
