#include "mol_descriptors.hpp"
#include "common.hpp"
#include "duckdb/common/assert.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/execution/expression_executor_state.hpp"
#include "duckdb/function/function_set.hpp"
#include "duckdb/main/extension_util.hpp"
#include "mol_formats.hpp"
#include "qed.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"

namespace duckdb_rdkit {

void mol_logp(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &binary_umbra_mol = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, float>(
      binary_umbra_mol, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto bmol = umbra_mol.GetBinaryMol();
        auto mol = rdkit_binary_mol_to_mol(bmol);
        double logp, _;
        RDKit::Descriptors::calcCrippenDescriptors(*mol, logp, _);
        return logp;
      });
}

void mol_qed(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &binary_umbra_mol = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, float>(
      binary_umbra_mol, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto bmol = umbra_mol.GetBinaryMol();
        auto mol = rdkit_binary_mol_to_mol(bmol);
        auto qed = QED();
        return qed.CalcQED(*mol);
      });
}

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

void mol_hbd(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &binary_umbra_mol = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, int32_t>(
      binary_umbra_mol, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto bmol = umbra_mol.GetBinaryMol();
        auto mol = rdkit_binary_mol_to_mol(bmol);
        return RDKit::Descriptors::calcNumHBD(*mol);
      });
}

void mol_hba(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &binary_umbra_mol = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, int32_t>(
      binary_umbra_mol, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto bmol = umbra_mol.GetBinaryMol();
        auto mol = rdkit_binary_mol_to_mol(bmol);
        return RDKit::Descriptors::calcNumHBA(*mol);
      });
}

void mol_num_rotatable_bonds(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &binary_umbra_mol = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, int32_t>(
      binary_umbra_mol, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto bmol = umbra_mol.GetBinaryMol();
        auto mol = rdkit_binary_mol_to_mol(bmol);
        return RDKit::Descriptors::calcNumRotatableBonds(*mol);
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

  ScalarFunctionSet set_mol_qed("mol_qed");
  set_mol_qed.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::FLOAT, mol_qed));
  ExtensionUtil::RegisterFunction(instance, set_mol_qed);

  ScalarFunctionSet set_mol_logp("mol_logp");
  set_mol_logp.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::FLOAT, mol_logp));
  ExtensionUtil::RegisterFunction(instance, set_mol_logp);

  ScalarFunctionSet set_mol_hbd("mol_hbd");
  set_mol_hbd.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::INTEGER, mol_hbd));
  ExtensionUtil::RegisterFunction(instance, set_mol_hbd);

  ScalarFunctionSet set_mol_hba("mol_hba");
  set_mol_hba.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::INTEGER, mol_hba));
  ExtensionUtil::RegisterFunction(instance, set_mol_hba);

  ScalarFunctionSet set_mol_num_rotatable_bonds("mol_num_rotatable_bonds");
  set_mol_num_rotatable_bonds.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::INTEGER, mol_num_rotatable_bonds));
  ExtensionUtil::RegisterFunction(instance, set_mol_num_rotatable_bonds);
}
} // namespace duckdb_rdkit
