#include "mol_descriptors.hpp"
#include "common.hpp"
#include "duckdb/common/assert.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "mol_formats.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"

namespace duckdb_rdkit {

void exactmw(DataChunk &args, ExpressionState &state, Vector &result) {
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

void RegisterDescriptorFunctions(DatabaseInstance &instance) {
  ScalarFunctionSet set_exactmw("exactmw");
  set_exactmw.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::FLOAT, exactmw));
  ExtensionUtil::RegisterFunction(instance, set_exactmw);
}
} // namespace duckdb_rdkit
