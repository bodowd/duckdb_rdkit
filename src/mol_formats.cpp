#include "mol_formats.hpp"
#include "common.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function_set.hpp"
#include "types.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <cstdint>
#include <memory>
#include <sys/types.h>

namespace duckdb_rdkit {
// Expects a SMILES string and returns a RDKit pickled molecule
std::unique_ptr<RDKit::ROMol> rdkit_mol_from_smiles(std::string s) {
  std::string smiles = s;
  std::unique_ptr<RDKit::ROMol> mol;
  try {
    mol.reset(RDKit::SmilesToMol(smiles));
  } catch (std::exception &e) {
    std::string msg = StringUtil::Format("%s", typeid(e).name());
    // TODO: throw a better exception. Right now, if it's not
    // a valid SMILES, it will break and then the whole db must be restarted
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

std::string serialize_umbra_mol(umbra_mol_t umbra_mol) {
  // encode a prefix and the binary mol into a single std::string
  // this can then be handed to the string_t constructor, and then the
  // constructor will take care of putting the prefix and the binary molecule
  // into the prefix and ptr fields in the `pointer` struct of string_t
  std::vector<char> buffer;

  buffer.insert(buffer.end(), reinterpret_cast<const char *>(&umbra_mol.prefix),
                reinterpret_cast<const char *>(&umbra_mol.prefix) +
                    sizeof(umbra_mol.prefix));
  buffer.insert(buffer.end(), umbra_mol.bmol.begin(), umbra_mol.bmol.end());

  return std::string(buffer.begin(), buffer.end());
}

std::string extract_bmol_from_umbra_mol(string_t buffer) {
  std::string bmol;

  // string_t::GetString() will get the data from the ptr to the string and
  // convert it to std::string
  auto prefix_size = string_t::PREFIX_BYTES;
  // string_t::GetString() will return the prefix and the bmol
  // extract just the bmol which is after 4 bytes of prefix
  // the total size of the string = prefix + bmol
  // so bmol_size = total size - prefix size
  auto bmol_size = buffer.GetSize() - string_t::PREFIX_BYTES;
  bmol.resize(bmol_size);
  std::memcpy(&bmol[0], &buffer.GetString()[prefix_size], bmol_size);

  std::cout << "extract_bmol_from_umbra_mol: " << std::endl;
  for (char b : bmol) {
    printf("%02x ", static_cast<unsigned char>(b));
  }

  return bmol;
}

uint32_t extract_prefix_from_umbra_mol(string_t buffer) {
  uint32_t prefix = Load<uint32_t>(const_data_ptr_cast(buffer.GetPrefix()));
  std::bitset<32> p(prefix);
  std::cout << p << '\n';
  return prefix;
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
        auto mol = rdkit_mol_from_smiles(smiles.GetString());
        auto pickled_mol = rdkit_mol_to_binary_mol(*mol);

        // add the meta data to the front of pickled mol and store the
        // buffer
        auto num_atoms = mol->getNumAtoms();
        auto num_bonds = mol->getNumBonds();
        auto amw = RDKit::Descriptors::calcAMW(*mol);
        auto num_rings = mol->getRingInfo()->numRings();
        std::cout << "First constructor call: " << std::endl;
        auto umbra_mol =
            umbra_mol_t(num_atoms, num_bonds, amw, num_rings, pickled_mol);
        std::cout << umbra_mol << std::endl;
        auto b_umbra_mol = serialize_umbra_mol(umbra_mol);
        std::cout << "b_umbra_mol: " << std::endl;
        for (auto i : b_umbra_mol) {
          printf("%02x ", static_cast<unsigned char>(i));
        }
        auto um = string_t(b_umbra_mol);
        std::cout << "\num.GetString()" << std::endl;
        for (char b : um.GetString()) {
          printf("%02x ", static_cast<unsigned char>(b));
        }

        return StringVector::AddString(result, um);
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
