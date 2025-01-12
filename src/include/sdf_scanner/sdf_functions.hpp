#pragma once
#include "duckdb/function/function_set.hpp"
#include "duckdb/function/replacement_scan.hpp"
namespace duckdb {

class SDFFunctions {
public:
  static vector<TableFunctionSet> GetTableFunctions();
  static unique_ptr<TableRef>
  ReadSDFReplacement(ClientContext &context, ReplacementScanInput &input,
                     optional_ptr<ReplacementScanData> data);

private:
  static TableFunctionSet GetReadSDFTableFunction();
  static TableFunctionSet GetReadSDFAutoTableFunction();
};
} // namespace duckdb
