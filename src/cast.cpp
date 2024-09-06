#include "common.hpp"
#include "duckdb/common/operator/cast_operators.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/function/cast/default_casts.hpp"
#include "duckdb/main/extension_util.hpp"
#include "mol_formats.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>

namespace duckdb_rdkit {

// This enables the user to insert into a Mol column by just writing the SMILES
// Duckdb will try to convert the string to a rdkit mol
// This is consistent with the RDKit Postgres cartridge behavior
void VarcharToMol(Vector &source, Vector &result, idx_t count) {
  UnaryExecutor::Execute<string_t, string_t>(
      source, result, count, [&](string_t smiles) {
        std::cout << "VarcharToMol called!" << std::endl;
        // this varchar is just a regular string, not a umbramol
        // Try to see if it is a SMILES
        auto mol = rdkit_mol_from_smiles(smiles.GetString());
        auto umbra_mol = get_umbra_mol_string(*mol);

        for (char b : umbra_mol) {
          printf("%02x ", static_cast<unsigned char>(b));
        }

        return StringVector::AddStringOrBlob(result, umbra_mol);
      });
}

bool VarcharToMolCast(Vector &source, Vector &result, idx_t count,
                      CastParameters &parameters) {
  VarcharToMol(source, result, count);
  return true;
}

void MolToVarchar(Vector &source, Vector &result, idx_t count) {
  UnaryExecutor::Execute<string_t, string_t>(
      source, result, count, [&](string_t b_umbra_mol) {
        // The input is a string_t coming from the duckdb internals.
        // The extension recognizes that this string_t is an
        // UmbraMol BLOB and will trigger this cast function.
        // Therefore, this function expects that the input
        // contains a string that has the format of umbra_mol_t.
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto bmol = umbra_mol.GetBinaryMol();

        auto rdkit_mol = rdkit_binary_mol_to_mol(bmol);
        auto smiles = rdkit_mol_to_smiles(*rdkit_mol);
        return StringVector::AddString(result, smiles);
      });
}

bool MolToVarcharCast(Vector &source, Vector &result, idx_t count,
                      CastParameters &parameters) {
  MolToVarchar(source, result, count);
  return true;
}

void RegisterCasts(DatabaseInstance &instance) {
  ExtensionUtil::RegisterCastFunction(instance, LogicalType::VARCHAR,
                                      ::duckdb_rdkit::Mol(),
                                      BoundCastInfo(VarcharToMolCast), 1);

  ExtensionUtil::RegisterCastFunction(instance, duckdb_rdkit::Mol(),
                                      LogicalType::VARCHAR,
                                      BoundCastInfo(MolToVarcharCast), 1);
}

} // namespace duckdb_rdkit
