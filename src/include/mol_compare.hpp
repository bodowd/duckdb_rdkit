#pragma once
#include "common.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "qed.hpp"
#include "umbra_mol.hpp"

namespace duckdb {
void RegisterCompareFunctions(ExtensionLoader &loader);

struct IsSubstructLocalState : public FunctionLocalState {
public:
  static unique_ptr<FunctionLocalState>
  Init(ExpressionState &state, const BoundFunctionExpression &expr,
       FunctionData *bind_data);
  void UpdateCache(string_t umbra_blob);
  unique_ptr<RDKit::ROMol> cached_mol;
};
} // namespace duckdb
