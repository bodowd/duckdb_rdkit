#include "mol_formats.hpp"
#include "common.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function_set.hpp"
#include "types.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Fingerprints/MACCS.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <memory>

namespace duckdb_rdkit {
// Expects a SMILES string and returns a RDKit pickled molecule
std::unique_ptr<RDKit::ROMol> rdkit_mol_from_smiles(std::string s) {
  std::string smiles = s;
  std::unique_ptr<RDKit::ROMol> mol;
  try {
    mol.reset(RDKit::SmilesToMol(smiles));
  } catch (std::exception &e) {
    std::string msg = StringUtil::Format("%s", typeid(e).name());
    // not sure if this is the right way to throw an error in duckdb
    throw FatalException(msg);
  }

  if (mol) {
    return mol;
  } else {
    std::string msg = StringUtil::Format("Could not convert %s to mol", smiles);
    throw FatalException(msg);
  }
}

// Serialize a molecule to binary using RDKit's MolPickler
std::string rdkit_mol_to_binary_mol(const RDKit::ROMol mol) {
  std::string buf;
  try {
    RDKit::MolPickler::pickleMol(mol, buf);
  } catch (...) {
    std::string msg = "Could not serialize mol to binary";
    throw FatalException(msg);
  }
  return buf;
}

// return std::string because it seems a little easier to work with with duckdb
// for example, StringVector. But perhaps there's a better way
std::string serialize_umbra_mol(umbra_mol_t umbra_mol) {
  std::vector<char> buffer;

  buffer.insert(buffer.end(),
                reinterpret_cast<const char *>(&umbra_mol.num_atoms),
                reinterpret_cast<const char *>(&umbra_mol.num_atoms) +
                    sizeof(umbra_mol.num_atoms));
  buffer.insert(buffer.end(),
                reinterpret_cast<const char *>(&umbra_mol.num_bonds),
                reinterpret_cast<const char *>(&umbra_mol.num_bonds) +
                    umbra_mol.NUM_BONDS_BYTES);
  buffer.insert(buffer.end(), reinterpret_cast<const char *>(&umbra_mol.amw),
                reinterpret_cast<const char *>(&umbra_mol.amw) +
                    umbra_mol.AMW_BYTES);
  buffer.insert(buffer.end(),
                reinterpret_cast<const char *>(&umbra_mol.num_rings),
                reinterpret_cast<const char *>(&umbra_mol.num_rings) +
                    umbra_mol.NUM_RINGS_BYTES);
  buffer.insert(buffer.end(),
                reinterpret_cast<const char *>(&umbra_mol.bmol_size),
                reinterpret_cast<const char *>(&umbra_mol.bmol_size) +
                    umbra_mol.BMOL_SIZE_BYTES);
  buffer.insert(buffer.end(), reinterpret_cast<const char *>(&umbra_mol.maccs),
                reinterpret_cast<const char *>(&umbra_mol.maccs) +
                    umbra_mol.MACCS_SIZE_BYTES);

  buffer.insert(buffer.end(), umbra_mol.bmol.begin(), umbra_mol.bmol.end());

  return std::string(buffer.begin(), buffer.end());
}

umbra_mol_t deserialize_umbra_mol(std::string buffer) {
  umbra_mol_t umbra_mol;
  size_t offset = 0;

  // Copy each member from the string buffer
  std::memcpy(&umbra_mol.num_atoms, &buffer[offset], umbra_mol.NUM_ATOMS_BYTES);
  offset += umbra_mol.NUM_ATOMS_BYTES;
  std::memcpy(&umbra_mol.num_bonds, &buffer[offset], umbra_mol.NUM_BONDS_BYTES);
  offset += umbra_mol.NUM_BONDS_BYTES;
  std::memcpy(&umbra_mol.amw, &buffer[offset], umbra_mol.AMW_BYTES);
  offset += umbra_mol.AMW_BYTES;
  std::memcpy(&umbra_mol.num_rings, &buffer[offset], umbra_mol.NUM_RINGS_BYTES);
  offset += umbra_mol.NUM_RINGS_BYTES;
  std::memcpy(&umbra_mol.bmol_size, &buffer[offset], umbra_mol.BMOL_SIZE_BYTES);
  offset += umbra_mol.BMOL_SIZE_BYTES;
  std::memcpy(&umbra_mol.maccs, &buffer[offset], umbra_mol.MACCS_SIZE_BYTES);
  offset += umbra_mol.MACCS_SIZE_BYTES;

  // std::vector<char> vec;
  // auto substring = buffer.substr(offset, offset + umbra_mol.bmol_size);
  // vec.insert(vec.end(), substring.begin(), substring.end());

  umbra_mol.bmol.resize(umbra_mol.bmol_size);
  // std::cout << "deserialize substr: " << std::endl;
  // for (auto i = offset; i < offset + umbra_mol.bmol_size; i++) {
  //   printf("%02x ", static_cast<unsigned char>(buffer[i]));
  // }
  std::memcpy(&umbra_mol.bmol[0], &buffer[offset], umbra_mol.bmol_size);

  return umbra_mol;
}

// Deserialize a binary mol to RDKit mol
std::unique_ptr<RDKit::ROMol> rdkit_binary_mol_to_mol(std::string bmol) {
  std::unique_ptr<RDKit::ROMol> mol(new RDKit::ROMol());
  RDKit::MolPickler::molFromPickle(bmol, *mol);

  return mol;
}

// Expects an RDKit pickled molecule and returns the SMILES of the molecule
std::string rdkit_mol_to_smiles(RDKit::ROMol mol) {
  std::string smiles = RDKit::MolToSmiles(mol);
  return smiles;
}

// An extension function callable from duckdb
// converts a serialized RDKit molecule to SMILES
//
// If there is a table mols with a column of type Mol
//
//
// select mol_to_smiles(*) from mols;
//
void mol_to_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &bmol = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, string_t>(
      bmol, result, count, [&](string_t bmol) {
        auto mol = rdkit_binary_mol_to_mol(bmol.GetString());
        auto smiles = rdkit_mol_to_smiles(*mol);
        return StringVector::AddString(result, smiles);
      });
}

void umbra_mol_from_smiles(DataChunk &args, ExpressionState &state,
                           Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles = args.data[0];
  auto count = args.size();

  UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
      smiles, result, count,
      [&](string_t smiles, ValidityMask &mask, idx_t idx) {
        try {
          auto mol = rdkit_mol_from_smiles(smiles.GetString());
          auto pickled_mol = rdkit_mol_to_binary_mol(*mol);

          // add the meta data to the front of pickled mol and store the
          // buffer
          auto num_atoms = mol->getNumAtoms();
          auto num_bonds = mol->getNumBonds();
          auto amw = RDKit::Descriptors::calcAMW(*mol);
          auto num_rings = mol->getRingInfo()->numRings();

          auto maccs_bit_vect = std::unique_ptr<ExplicitBitVect>{
              RDKit::MACCSFingerprints::getFingerprintAsBitVect(*mol)};

          auto umbra_mol = umbra_mol_t(num_atoms, num_bonds, amw, num_rings,
                                       std::move(maccs_bit_vect), pickled_mol);

          auto b_umbra_mol = serialize_umbra_mol(umbra_mol);

          return StringVector::AddString(result, b_umbra_mol);
        } catch (...) {
          mask.SetInvalid(idx);
          return string_t();
        }
      });
}

// An extension function callable from duckdb
// converts a SMILES all the way to the serialized version of the RDKit mol
// returns NULL if conversion fails
//
// select mol_from_smiles('CC');
//
void mol_from_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles = args.data[0];
  auto count = args.size();

  UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
      smiles, result, count,
      [&](string_t smiles, ValidityMask &mask, idx_t idx) {
        try {
          auto mol = rdkit_mol_from_smiles(smiles.GetString());
          auto pickled_mol = rdkit_mol_to_binary_mol(*mol);
          return StringVector::AddString(result, pickled_mol);
        } catch (...) {
          mask.SetInvalid(idx);
          return string_t();
        }
      });
}

void RegisterFormatFunctions(DatabaseInstance &instance) {
  // Register scalar functions
  ScalarFunctionSet mol_from_smiles_set("mol_from_smiles");
  mol_from_smiles_set.AddFunction(ScalarFunction(
      {LogicalType::VARCHAR}, duckdb_rdkit::Mol(), mol_from_smiles));
  ExtensionUtil::RegisterFunction(instance, mol_from_smiles_set);

  ScalarFunctionSet mol_to_smiles_set("mol_to_smiles");
  mol_to_smiles_set.AddFunction(ScalarFunction(
      {duckdb_rdkit::Mol()}, LogicalType::VARCHAR, mol_to_smiles));
  ExtensionUtil::RegisterFunction(instance, mol_to_smiles_set);

  ScalarFunctionSet umbra_mol_from_smiles_set("umbra_mol_from_smiles");
  umbra_mol_from_smiles_set.AddFunction(ScalarFunction(
      {LogicalType::VARCHAR}, duckdb_rdkit::UmbraMol(), umbra_mol_from_smiles));
  ExtensionUtil::RegisterFunction(instance, umbra_mol_from_smiles_set);
}

} // namespace duckdb_rdkit
