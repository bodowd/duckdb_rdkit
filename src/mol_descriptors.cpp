#include "mol_descriptors.hpp"
#include "common.hpp"
#include "duckdb/common/assert.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/execution/expression_executor_state.hpp"
#include "mol_formats.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"

namespace duckdb_rdkit {

void mol_amw(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &binary_umbra_mol = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, float>(
      binary_umbra_mol, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto bmol = umbra_mol.GetBinaryMol();
        auto mol = rdkit_binary_mol_to_mol(bmol);
        return RDKit::Descriptors::calcAMW(*mol);
      });
}

void mol_exactmw(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &binary_umbra_mol = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, float>(
      binary_umbra_mol, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto bmol = umbra_mol.GetBinaryMol();
        auto mol = rdkit_binary_mol_to_mol(bmol);
        return RDKit::Descriptors::calcExactMW(*mol);
      });
}

void mol_tpsa(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &binary_umbra_mol = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, float>(
      binary_umbra_mol, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto bmol = umbra_mol.GetBinaryMol();
        auto mol = rdkit_binary_mol_to_mol(bmol);
        return RDKit::Descriptors::calcTPSA(*mol);
      });
}

void RegisterDescriptorFunctions(DatabaseInstance &instance) {
  ScalarFunctionSet set_mol_amw("mol_amw");
  set_mol_amw.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::FLOAT, mol_amw));
  ExtensionUtil::RegisterFunction(instance, set_mol_amw);

  ScalarFunctionSet set_mol_exactmw("mol_exactmw");
  set_mol_exactmw.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::FLOAT, mol_exactmw));
  ExtensionUtil::RegisterFunction(instance, set_mol_exactmw);

  ScalarFunctionSet set_mol_tpsa("mol_tpsa");
  set_mol_tpsa.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::FLOAT, mol_tpsa));
  ExtensionUtil::RegisterFunction(instance, set_mol_tpsa);
}
} // namespace duckdb_rdkit
