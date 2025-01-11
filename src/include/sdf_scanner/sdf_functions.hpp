#pragma once
#include "duckdb/function/function_set.hpp"
namespace duckdb {

class SDFFunctions {
public:
  static vector<TableFunctionSet> GetTableFunctions();

private:
  static TableFunctionSet GetReadSDFTableFunction();
};
} // namespace duckdb
