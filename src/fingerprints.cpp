#include "fingerprints.hpp"
#include "common.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "mol_formats.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"
#include <GraphMol/Fingerprints/MorganFingerprints.h>

namespace duckdb_rdkit {

std::string get_morganbv_fp(const RDKit::ROMol &mol, int radius) {
  std::cout << "get_morganbv_fp" << std::endl;
  std::vector<std::uint32_t> invars(mol.getNumAtoms());
  // try {
  RDKit::MorganFingerprints::getConnectivityInvariants(mol, invars, true);
  auto res = std::unique_ptr<ExplicitBitVect>(
      RDKit::MorganFingerprints::getFingerprintAsBitVect(mol, radius, 1024,
                                                         &invars));
  if (res) {
    std::cout << "RES!" << std::endl;
    for (char b : res->toString()) {
      printf("%02x ", static_cast<unsigned char>(b));
    }
    std::cout << "" << std::endl;
    for (auto i = 0; i < res->size(); i++) {
      printf("%d", res->getBit(i));
    }
  }
  // } catch (std::exception &e) {
  //   std::string msg =
  //       StringUtil::Format("error in get_morganbv_fp: %s", typeid(e).name());
  //   throw InvalidInputException(msg);
  // }
  std::cout << "returning ... " << std::endl;
  return res->toString();
}

void morganbv_fp(DataChunk &args, ExpressionState &state, Vector &result) {
  // this expects a Mol as the argument
  // `umbra_mol` is a special mol object in duckdb_rdkit. It has more
  // information in addition to the rdkit mol (see umbra_mol.hpp for details)
  auto &umbra_mol = args.data[0];
  std::cout << "in morganbv_fp: " << std::endl;
  auto count = args.size();

  UnaryExecutor::Execute<string_t, string_t>(
      umbra_mol, result, count, [&](umbra_mol_t umbra_mol) {
        // std::cout << "umbra_mol" << std::endl;
        // for (char b : umbra_mol.GetString()) {
        //   printf("%02x ", static_cast<unsigned char>(b));
        // }
        // std::cout << "" << std::endl;

        auto bmol = umbra_mol.GetBinaryMol();
        auto mol = rdkit_binary_mol_to_mol(bmol);
        auto morganbv_fp = get_morganbv_fp(*mol, 2);
        // std::cout << "umbra_mol bmol: " << std::endl;
        for (char b : morganbv_fp) {
          printf("%02x ", static_cast<unsigned char>(b));
        }
        return StringVector::AddStringOrBlob(result, morganbv_fp);
      });
}

void RegisterFingerprintFunctions(DatabaseInstance &instance) {
  ScalarFunctionSet morganbv_fp_set("morganbv_fp");
  morganbv_fp_set.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::VARCHAR, morganbv_fp));

  ExtensionUtil::RegisterFunction(instance, morganbv_fp_set);
}

} // namespace duckdb_rdkit
