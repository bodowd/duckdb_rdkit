#include "qed.hpp"
#include "mol_formats.hpp"
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace duckdb_rdkit {

std::vector<RDKit::RWMol> smarts2mols(std::vector<std::string> smarts) {
  std::vector<RDKit::RWMol> mols;
  for (auto s : smarts) {
    std::unique_ptr<RDKit::RWMol> mol;
    mol.reset(RDKit::SmartsToMol(s));
    mols.push_back(*mol);
  }
  return mols;
}

double QED::CalcADS(float x, std::string adsParameterKey) {
  auto p = QED::adsParameters[adsParameterKey];
  auto exp1 = 1 + std::exp(-1 * (x - p.C + p.D / 2) / p.E);
  auto exp2 = 1 + std::exp(-1 * (x - p.C - p.D / 2) / p.F);
  auto dx = p.A + p.B / exp1 * (1 - 1 / exp2);
  return dx / p.DMAX;
}

QED::QEDproperties
QED::CalcProperties(const RDKit::ROMol &mol,
                    std::vector<RDKit::RWMol> acceptorMols,
                    std::unique_ptr<RDKit::RWMol> aliphaticRingMol,
                    std::vector<RDKit::RWMol> alertMols) {
  auto molWithoutHs = RDKit::MolOps::removeHs(mol);
  auto amw = RDKit::Descriptors::calcAMW(mol);
  double crippenLogP = 0;
  double _mr = 0;
  RDKit::Descriptors::calcCrippenDescriptors(mol, crippenLogP, _mr);

  // find all hydrogen bond acceptors
  RDKit::MatchVectType matchVect;
  auto hba = 0;
  for (auto m : acceptorMols) {
    bool match = RDKit::SubstructMatch(mol, m, matchVect);
    hba += matchVect.size();
  }

  auto hbd = RDKit::Descriptors::calcNumHBD(mol);
  auto psa = RDKit::Descriptors::calcTPSA(mol);
  auto rotb = RDKit::Descriptors::calcNumRotatableBonds(mol, true);
  auto withoutAliphaticRings = RDKit::deleteSubstructs(mol, *aliphaticRingMol);
  auto arom = RDKit::MolOps::findSSSR(*withoutAliphaticRings);

  auto alerts = 0;
  for (auto &m : alertMols) {
    bool match = RDKit::SubstructMatch(mol, m, matchVect);
    if (match) {
      alerts += 1;
    }
  }

  return QEDproperties(amw, crippenLogP, hba, hbd, psa, rotb, arom, alerts);
}

float QED::CalcQED(const RDKit::ROMol &mol) {
  auto properties = CalcProperties(mol, acceptorMols,
                                   std::move(aliphaticRingsMol), alertMols);
  float sumOfWeightedADSValues = 0.0;
  float sumOfWeights = 0.0;

  for (const auto &[k, v] : properties.data) {
    sumOfWeightedADSValues +=
        (WEIGHT_MEAN.data[k] * std::log(CalcADS(properties.data[k], k)));
    sumOfWeights += WEIGHT_MEAN.data[k];
  }

  return std::exp(sumOfWeightedADSValues / sumOfWeights);
}

} // namespace duckdb_rdkit
