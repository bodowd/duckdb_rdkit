#include "common.hpp"
#include "duckdb/common/string_util.hpp"
#include "types.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

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

          // add the meta data to the front of pickled mol and store the buffer
          auto num_atoms = mol->getNumAtoms();
          auto num_bonds = mol->getNumBonds();
          auto amw = RDKit::Descriptors::calcAMW(*mol);
          auto num_rings = mol->getRingInfo()->numRings();

          auto pickled_mol = rdkit_mol_to_binary_mol(*mol);

          auto umbra_mol =
              UmbraMol(num_atoms, num_bonds, amw, num_rings, pickled_mol);

          std::cout << "\numbra mol" << std::endl;
          std::cout << "num_atoms: " << umbra_mol.num_atoms << std::endl;
          std::cout << "num_bonds: " << umbra_mol.num_bonds << std::endl;
          std::cout << "amw: " << umbra_mol.amw << std::endl;
          std::cout << "num_rings: " << umbra_mol.num_rings << std::endl;
          std::cout << "bmol_size: " << umbra_mol.bmol_size << std::endl;
          // std::cout << "binary_mol: " << umbra_mol.bmol_ptr.get() <<
          // std::endl;
          std::vector<char> pbmol;
          for (size_t i = 0; i < umbra_mol.bmol_size; ++i) {
            auto bmol = umbra_mol.bmol;
            printf("%02x ", static_cast<unsigned char>(bmol[i]));
            pbmol.push_back(bmol[i]);
          }
          printf("\n");

          std::cout << "vector: " << std::endl;
          for (auto i : pbmol) {
            printf("%02x ", static_cast<unsigned char>(i));
          }

          std::cout << "pickled mol: " << std::endl;
          for (auto i : pickled_mol) {
            printf("%02x ", static_cast<unsigned char>(i));
          }

          // auto serialized = umbra_mol.serialize();
          // std::cout << "serialized" << std::endl;
          // for (auto i : serialized) {
          //   // printf("%02x", i);
          //   printf("%02x ", static_cast<unsigned char>(i));
          // }

          // NOTE: it seems that RDKit expects a std::string for the depickling
          // When I change the type of the binary mol that i keep in umbra_mol
          // it works. Perhaps RDKit uses some metadata from the std::string
          // that is missing when just getting the .data() from the std::string?

          // TODO: convert the whole serialized thing into a string so that it
          // can be passed to a StringVector
          return StringVector::AddString(result, umbra_mol.bmol);
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
}

} // namespace duckdb_rdkit
