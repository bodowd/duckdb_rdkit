#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/execution/expression_executor_state.hpp"
#include "duckdb/function/function_set.hpp"
#include <exception>
#define DUCKDB_EXTENSION_MAIN

#include "duckdb/common/exception.hpp"
#include "duckdb/common/operator/cast_operators.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension_util.hpp"
#include "duckdb_rdkit_extension.hpp"
#include "types.hpp"
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>
#include <string>

namespace duckdb {

// Serialize a molecule to binary using RDKit's MolPickler
//
string_t rdkit_mol_to_binary_mol(const RDKit::ROMol &mol) {
  std::string buf;
  try {
    RDKit::MolPickler::pickleMol(mol, buf, RDKit::PicklerOps::AllProps);
  } catch (...) {
    std::string msg = "Could not serialize mol to binary";
    throw Exception(msg);
  }
  return buf;
}

RDKit::ROMol rdkit_binary_mol_to_mol(std::string &bmol) {
  RDKit::ROMol mol;
  RDKit::MolPickler::molFromPickle(bmol, mol, RDKit::PicklerOps::AllProps);
  return mol;
}

static std::string rdkit_mol_from_smiles(std::string s) {
  /* build the molecule blob repr from a text string */
  // std::string smiles = s;
  // RDKit::ROMol *mol;
  // try {
  //   mol = RDKit::SmilesToMol(smiles);
  // } catch (std::exception &e) {
  //   // std::string msg = StringUtil::Format("Could not convert %s to mol",
  //   // smiles);
  //   std::string msg = StringUtil::Format("%s", typeid(e).name());
  //   // not sure if this is the right way to throw an error in duckdb
  //   throw Exception(msg);
  // }
  //
  // if (mol) {
  //   // serialize the mol
  //   // std::string bmol = rdkit_mol_to_binary_mol(*mol);
  //   std::string buf;
  //   RDKit::MolPickler::pickleMol(*mol, buf);
  //   return buf;
  // } else {
  //   std::string msg = StringUtil::Format("Could not convert %s to mol",
  //   smiles); throw Exception(msg);
  // }
  //
  RDKit::ROMol *mol = RDKit::SmilesToMol(s);
  std::string pickle;
  RDKit::MolPickler::pickleMol(*mol, pickle);
  std::cout << pickle << std::endl;
  return pickle;
}

static std::string rdkit_mol_to_smiles(std::string bmol) {
  // RDKit::ROMol mol = rdkit_binary_mol_to_mol(bmol);
  std::unique_ptr<RDKit::ROMol> mol(new RDKit::ROMol());
  RDKit::MolPickler::molFromPickle(bmol, *mol);
  std::string smiles = RDKit::MolToSmiles(*mol);
  return smiles;
}

inline void DuckdbRdkitScalarFun(DataChunk &args, ExpressionState &state,
                                 Vector &result) {
  auto &smiles_vector = args.data[0];
  UnaryExecutor::Execute<string_t, string_t>(
      smiles_vector, result, args.size(), [&](string_t s) {
        // auto b = rdkit_mol_from_smiles(s.GetString());
        RDKit::ROMol *mol = RDKit::SmilesToMol(s.GetString());
        std::string pickle;
        RDKit::MolPickler::pickleMol(*mol, pickle);
        std::string msg = StringUtil::Format("%s", pickle);

        HandleCastError::AssignError(msg, &msg);
        return StringVector::AddString(result, "DuckdbRdkit " + msg + "🐥");
        ;
      });
}

void mol_from_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, string_t>(
      smiles, result, count, [&](string_t smiles) {
        auto binary_mol = rdkit_mol_from_smiles(smiles.GetString());
        // TODO: seems to convert this to rdkit canonical
        // but using mol_to_smiles(mol_from_smiles('CC')) in duckdb doesn't
        // work
        // somehow the data is not in a correct binary format anymore
        // maybe something needs to be cast?
        // I think I'm casint to loicaltype blob called 'Mol' in duckdb
        // and i need to cast that back to a string before it can be
        // deserialized?
        // return rdkit_mol_to_smiles(binary_mol);
        // return binary_mol;
        RDKit::ROMol *mol = RDKit::SmilesToMol(smiles.GetString());
        std::string pickle;
        RDKit::MolPickler::pickleMol(*mol, pickle);
        return StringVector::AddString(result, pickle);
      });
}

void mol_to_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, string_t>(
      smiles, result, count, [&](string_t bmol) {
        auto smiles = rdkit_mol_to_smiles(bmol.GetString());
        return StringVector::AddString(result, smiles);
        // RDKit::ROMol mol;
        // RDKit::MolPickler::molFromPickle(bmol.GetString(),
        // mol,
        //                                  RDKit::PicklerOps::AllProps);
        // return mol;
        //

        // auto smol = rdkit_mol_from_smiles(bmol.GetString());
        // auto smol = RDKit::SmilesToMol(bmol.GetString());
        // return RDKit::MolToSmiles(*smol);
        // RDKit::ROMol mol;
        // RDKit::MolPickler::molFromPickle(smol, mol,
        //                                  RDKit::PicklerOps::AllProps);

        // return
        // rdkit_mol_from_smiles(smol);
      });
}

// void MolToVarchar(Vector &source, Vector &result, idx_t count) {
//   using VARCHAR_TYPE = PrimitiveType<string_t>;
//   using MOL_TYPE = StructTypeUnary<long>;
//   GenericExecutor::ExecuteUnary<MOL_TYPE, VARCHAR_TYPE>(
//       source, result, count, [&](MOL_TYPE &mol) {
//         return StringVector::AddString(result,
//                                        StringUtil::Format("Mol (%s)", mol));
//       });
// }
//
// static bool MolToVarcharCast(Vector &source, Vector &result, idx_t count,
//                              CastParameters &parameters) {
//   MolToVarchar(source, result, count);
//   return true;
// }

// bool CastVarcharToMol(Vector &source, Vector &result, idx_t count,
//                       CastParameters &parameters) {
//   auto &entries = StructVector::GetEntries(source);
//   auto &smiles = entries[0];
//
//   auto &result_children = StructVector::GetEntries(result);
//   auto &mol = result_children[0];
//
//   if (count == 1) {
//     result.SetVectorType(VectorType::CONSTANT_VECTOR);
//   }
//   return true;
// }

// bool CastMolToVarchar(Vector &source, Vector &result, idx_t count,
//                       CastParameters &parameters) {
//   auto &entries = StructVector::GetEntries(source);
//   auto &mol = entries[0];
//
//   auto &result_children = StructVector::GetEntries(result);
//   auto &smiles = result_children[0];
//
//   if (count == 1) {
//     result.SetVectorType(VectorType::CONSTANT_VECTOR);
//   }
//
//   UnaryExecutor::Execute<string_t, string_t>(
//       source, result, count,
//       [&](string_t bmol) { return rdkit_mol_to_smiles(bmol.GetString()); });
//   return true;
// }

static void LoadInternal(DatabaseInstance &instance) {
  // Register Mol type
  ExtensionUtil::RegisterType(instance, "Mol", duckdb_rdkit::Mol());
  // Register a scalar function
  auto duckdb_rdkit_scalar_function =
      ScalarFunction("duckdb_rdkit", {LogicalType::VARCHAR},
                     LogicalType::VARCHAR, DuckdbRdkitScalarFun);
  ExtensionUtil::RegisterFunction(instance, duckdb_rdkit_scalar_function);

  ScalarFunctionSet mol_from_smiles_set("mol_from_smiles");
  mol_from_smiles_set.AddFunction(ScalarFunction(
      {LogicalType::VARCHAR}, duckdb_rdkit::Mol(), mol_from_smiles));
  ExtensionUtil::RegisterFunction(instance, mol_from_smiles_set);

  ScalarFunctionSet mol_to_smiles_set("mol_to_smiles");
  mol_to_smiles_set.AddFunction(ScalarFunction(
      {duckdb_rdkit::Mol()}, LogicalType::VARCHAR, mol_to_smiles));
  ExtensionUtil::RegisterFunction(instance, mol_to_smiles_set);

  // ExtensionUtil::RegisterCastFunction(instance, duckdb_rdkit::Mol(),
  //                                     LogicalType::VARCHAR,
  //                                     DefaultCasts::ReinterpretCast, 1);
  // ExtensionUtil::RegisterCastFunction(instance, LogicalType::VARCHAR,
  //                                     duckdb_rdkit::Mol(), CastVarcharToMol);

  // ExtensionUtil::RegisterCastFunction(instance, duckdb_rdkit::Mol(),
  //                                     LogicalType::VARCHAR,
  //                                     CastMolToVarchar);
}

void DuckdbRdkitExtension::Load(DuckDB &db) { LoadInternal(*db.instance); }

std::string DuckdbRdkitExtension::Name() { return "duckdb_rdkit"; }

} // namespace duckdb

extern "C" {

DUCKDB_EXTENSION_API void duckdb_rdkit_init(duckdb::DatabaseInstance &db) {
  duckdb::DuckDB db_wrapper(db);
  db_wrapper.LoadExtension<duckdb::DuckdbRdkitExtension>();
}

DUCKDB_EXTENSION_API const char *duckdb_rdkit_version() {
  return duckdb::DuckDB::LibraryVersion();
}
}

#ifndef DUCKDB_EXTENSION_MAIN
#error DUCKDB_EXTENSION_MAIN not defined
#endif
