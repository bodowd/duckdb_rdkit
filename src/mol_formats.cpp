#include "mol_formats.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/helper.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/unique_ptr.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "types.hpp"
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <memory>

namespace duckdb_rdkit {
// Expects a SMILES string and returns a RDKit pickled molecule
std::unique_ptr<RDKit::ROMol> rdkit_mol_from_smiles(std::string smiles) {
  std::cout << "SMILES in rdkit_mol_from_smiles: " << typeid(smiles).name()
            << std::endl;
  std::unique_ptr<RDKit::ROMol> mol;
  try {
    std::cout << RDKit::SmilesToMol(smiles) << std::endl;
    mol.reset(RDKit::SmilesToMol(smiles));
  } catch (std::exception &e) {
    std::string msg = StringUtil::Format("%s", typeid(e).name());
    // not sure if this is the right way to throw an error in duckdb
    throw FatalException(msg);
  }

  if (mol) {
    return std::move(mol);
  } else {
    string msg = StringUtil::Format("Could not convert %s to mol", smiles);
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
        // TODO: convert mol_to* and mol_from* to return bool so that
        // a mask can be set if not valid
        // like so:
        //
        // HandleCastError::AssignError(msg, parameters);
        // mask.SetInvalid(idx);
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
  auto &func_expr = state.expr.Cast<BoundFunctionExpression>();
  auto &info = func_expr.bind_info->Cast<BindData>();

  auto bound_smiles = info.smiles;
  std::cout << "bound smiles: " << RDKit::SmilesToMol(bound_smiles)
            << std::endl;

  UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
      smiles, result, count,
      [&](string_t smiles, ValidityMask &mask, idx_t idx) {
        try {
          // auto mol = rdkit_mol_from_smiles(smiles.GetString());
          auto mol = rdkit_mol_from_smiles(bound_smiles);
          auto pickled_mol = rdkit_mol_to_binary_mol(*mol);
          return StringVector::AddString(result, pickled_mol);
        } catch (...) {
          mask.SetInvalid(idx);
          return string_t();
        }
      });
}

BindData::BindData(const BindData &other) : smiles(other.smiles) {}

BindData::BindData(const string &smiles) : smiles(smiles) {}

bool BindData::Equals(const FunctionData &other_p) const {
  auto &other = other_p.Cast<const BindData>();
  return other.smiles == smiles;
  // return mol_cmp(*mol, *other.mol);
}

unique_ptr<FunctionData> BindData::Copy() const {
  return make_uniq<BindData>(smiles);
}

unique_ptr<FunctionData> Bind(ClientContext &context,
                              ScalarFunction &bound_function,
                              vector<unique_ptr<Expression>> &arguments) {
  D_ASSERT(bound_function.arguments.size() == 1);
  std::cout << "arguments[0]" << std::endl;
  auto smiles = arguments[0]->ToString();
  std::cout << arguments[0]->ToString() << std::endl;
  // auto bind_data = make_uniq<BindData>(smiles);
  // bind_data->mol = std::move(rdkit_mol_from_smiles(smiles));

  // if (arguments[0]->IsFoldable()) {
  //   const auto path_val =
  //       ExpressionExecutor::EvaluateScalar(context, *arguments[0]);
  // }
  // std::cout << "Arguments in Bind: " << smiles << std::endl;
  // // const auto path_val =
  // //     ExpressionExecutor::EvaluateScalar(context, *arguments[0]);
  // // std::cout << "PATH VAL in BIND: " << path_val << std::endl;
  // bound_function.arguments[0] = LogicalType::VARCHAR;
  // auto molfromsmiles = rdkit_mol_from_smiles(smiles.GetString());
  // std::cout << "# of atoms in mol in Bind: " << molfromsmiles->getNumAtoms()
  //           << std::endl;
  //
  // auto mol = make_uniq<RDKit::ROMol>(rdkit_mol_from_smiles(smiles));
  return make_uniq<BindData>(smiles);
}

void RegisterFormatFunctions(DatabaseInstance &instance) {
  // Register scalar functions
  ScalarFunctionSet mol_from_smiles_set("mol_from_smiles");
  mol_from_smiles_set.AddFunction(ScalarFunction(
      {LogicalType::VARCHAR}, duckdb_rdkit::Mol(), mol_from_smiles, Bind));
  ExtensionUtil::RegisterFunction(instance, mol_from_smiles_set);

  ScalarFunctionSet mol_to_smiles_set("mol_to_smiles");
  mol_to_smiles_set.AddFunction(ScalarFunction(
      {duckdb_rdkit::Mol()}, LogicalType::VARCHAR, mol_to_smiles));
  ExtensionUtil::RegisterFunction(instance, mol_to_smiles_set);
}

} // namespace duckdb_rdkit
