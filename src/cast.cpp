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
        auto mol = rdkit_mol_from_smiles(smiles.GetString());
        auto pickled_mol = rdkit_mol_to_binary_mol(*mol);
        return StringVector::AddString(result, pickled_mol);
      });
}

bool VarcharToMolCast(Vector &source, Vector &result, idx_t count,
                      CastParameters &parameters) {
  VarcharToMol(source, result, count);
  return true;
}

void MolToVarchar(Vector &source, Vector &result, idx_t count) {
  UnaryExecutor::Execute<string_t, string_t>(
      source, result, count, [&](string_t bmol) {
        auto mol = rdkit_binary_mol_to_mol(bmol.GetString());
        auto smiles = rdkit_mol_to_smiles(*mol);
        return StringVector::AddString(result, smiles);
      });
}

bool MolToVarcharCast(Vector &source, Vector &result, idx_t count,
                      CastParameters &parameters) {
  MolToVarchar(source, result, count);
  return true;
}

void UmbraMolToVarchar(Vector &source, Vector &result, idx_t count) {
  UnaryExecutor::Execute<string_t, string_t>(
      source, result, count, [&](string_t b_umbra_mol) {
        // The input is a string_t coming from the duckdb system.
        // The extension recognizes that this string_t is an
        // UmbraMol BLOB and will trigger this cast function.
        // Therefore, this function expects that the input
        // contains a string that has the format of umbra_mol_t

        // std::cout << "\nUmbraMolToVarchar" << std::endl;
        // for (char b : b_umbra_mol.GetString()) {
        //   printf("%02x ", static_cast<unsigned char>(b));
        // }

        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto bmol = umbra_mol.GetBinaryMol();
        // std::cout << "\nbmol from umbramol to varchar" << std::endl;
        // for (char b : bmol) {
        //   printf("%02x ", static_cast<unsigned char>(b));
        // }
        // std::cout << "UMBRA MOL GET STRING: " << std::endl;
        // for (char b : umbra_mol.GetString()) {
        //   printf("%02x ", static_cast<unsigned char>(b));
        // }

        auto rdkit_mol = rdkit_binary_mol_to_mol(bmol);
        auto smiles = rdkit_mol_to_smiles(*rdkit_mol);
        return StringVector::AddString(result, smiles);
      });
}

bool UmbraMolToVarcharCast(Vector &source, Vector &result, idx_t count,
                           CastParameters &parameters) {
  UmbraMolToVarchar(source, result, count);
  return true;
}

void VarcharToUmbraMol(Vector &source, Vector &result, idx_t count) {
  UnaryExecutor::Execute<string_t, string_t>(
      source, result, count, [&](string_t smiles) {
        std::cout << "VARCHAR to UMBRA MOL" << std::endl;
        // this varchar is just a regular string, not a umbramol
        auto mol = rdkit_mol_from_smiles(smiles.GetString());
        // add the meta data to the front of pickled mol and store the buffer
        auto num_atoms = mol->getNumAtoms();
        auto num_bonds = mol->getNumBonds();
        auto amw = RDKit::Descriptors::calcAMW(*mol);
        auto num_rings = mol->getRingInfo()->numRings();

        auto pickled_mol = rdkit_mol_to_binary_mol(*mol);
        auto umbra_mol = get_umbra_mol_string(num_atoms, num_bonds, amw,
                                              num_rings, pickled_mol);

        return StringVector::AddStringOrBlob(result, umbra_mol);
      });
}

bool VarcharToUmbraMolCast(Vector &source, Vector &result, idx_t count,
                           CastParameters &parameters) {
  VarcharToUmbraMol(source, result, count);
  return true;
}

void RegisterCasts(DatabaseInstance &instance) {
  ExtensionUtil::RegisterCastFunction(instance, LogicalType::VARCHAR,
                                      ::duckdb_rdkit::Mol(),
                                      BoundCastInfo(VarcharToMolCast), 1);

  ExtensionUtil::RegisterCastFunction(instance, duckdb_rdkit::Mol(),
                                      LogicalType::VARCHAR,
                                      BoundCastInfo(MolToVarcharCast), 1);

  ExtensionUtil::RegisterCastFunction(instance, duckdb_rdkit::UmbraMol(),
                                      LogicalType::VARCHAR,
                                      BoundCastInfo(UmbraMolToVarcharCast), 1);

  ExtensionUtil::RegisterCastFunction(instance, LogicalType::VARCHAR,
                                      duckdb_rdkit::UmbraMol(),
                                      BoundCastInfo(VarcharToUmbraMolCast), 1);
}

} // namespace duckdb_rdkit
