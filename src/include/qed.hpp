#pragma once

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <memory>
#include <unordered_map>

namespace duckdb_rdkit {

// Quantitative estimation of drug-likeness (QED) is implemented according to
// the RDKit python implementation.
//
// QED paper: Bickerton, G.R.; Paolini, G.V.; Besnard, J.; Muresan, S.; Hopkins,
// A.L. (2012) 'Quantifying the chemical beauty of drugs', Nature Chemistry, 4,
// 90-98
class QED {
public:
  QED() {
    acceptorMols = smarts2mols(acceptorSmarts);
    alertMols = smarts2mols(structuralAlertSmarts);
    aliphaticRingsMol.reset(RDKit::SmartsToMol(aliphaticRingSmarts));
  }

  float CalcQED(const RDKit::ROMol &mol);

private:
  struct ADSparameter {
    double A, B, C, D, E, F, DMAX;
  };

  struct QEDproperties {
    std::unordered_map<std::string, float> data;

    // This constructor populates the internal map that stores the calculated
    // properties for the molecule we are calculating the QED for.
    // these fields need to be floats because the weights can be floats
    QEDproperties(float mw, float alogp, float hba, float hbd, float psa,
                  float rotb, float arom, float alerts) {
      // molecular weight
      data["MW"] = mw;
      // Crippen logp
      data["ALOGP"] = alogp;
      // number of hydrogen bond acceptors
      data["HBA"] = hba;
      // number of hydrogen bond donors
      data["HBD"] = hbd;
      // polar surface area
      data["PSA"] = psa;
      // number of rotatable bonds
      data["ROTB"] = rotb;
      // number of aromatic rings
      data["AROM"] = arom;
      // alerts = from the QED paper:"a substructure search was performed
      // against each drug using a curated reference set of 94 functional
      // moieties that are potentially mutagenic, reactive or have unfavourable
      // pharmacokinetic properties"
      data["ALERTS"] = alerts;
    };
  };

  // Helper function to convert a vector of SMARTS to a vector of RDKit
  // molecules
  std::vector<RDKit::RWMol> smarts2mols(std::vector<std::string> smarts);
  //  Asymmetric Double Sigmoidal (ADS) function parameters used to model the
  //  histogram. Parameters provided by the paper (see top of file)
  std::unordered_map<std::string, ADSparameter> adsParameters = {
      {"MW",
       {2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677,
        65.37051707, 104.9805561}},
      {"ALOGP",
       {3.172690585, 137.8624751, 2.534937431, 4.581497897, 0.822739154,
        0.576295591, 131.3186604}},
      {"HBA",
       {2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953,
        1.300669958, 148.7763046}},
      {"HBD",
       {1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843,
        0.920922555, 258.1632616}},
      {"PSA",
       {1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824,
        28.51324732, 104.5686167}},
      {"ROTB",
       {0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166,
        2.758063707, 105.4420403}},
      {"AROM",
       {3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384,
        0.375760881, 312.3372610}},
      {"ALERTS",
       {0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477,
        0.875193782, 417.7253140}}};
  std::vector<RDKit::RWMol> acceptorMols;
  std::vector<RDKit::RWMol> alertMols;
  std::unique_ptr<RDKit::RWMol> aliphaticRingsMol;
  std::string aliphaticRingSmarts = "[$([A;R][!a])]";
  std::vector<std::string> acceptorSmarts = {"[oH0;X2]",
                                             "[OH1;X2;v2]",
                                             "[OH0;X2;v2]",
                                             "[OH0;X1;v2]",
                                             "[O-;X1]",
                                             "[SH0;X2;v2]",
                                             "[SH0;X1;v2]",
                                             "[S-;X1]",
                                             "[nH0;X2]",
                                             "[NH0;X1;v3]",
                                             "[$([N;+0;X3;v3]);!$(N[C,S]=O)]"};
  // clang-format off
  std::vector<std::string> structuralAlertSmarts = {
 "*1[O,S,N]*1", "[S,C](=[O,S])[F,Br,Cl,I]", "[CX4][Cl,Br,I]", "[#6]S(=O)(=O)O[#6]",
  "[$([CH]),$(CC)]#CC(=O)[#6]", "[$([CH]),$(CC)]#CC(=O)O[#6]", "n[OH]",
  "[$([CH]),$(CC)]#CS(=O)(=O)[#6]", "C=C(C=O)C=O", "n1c([F,Cl,Br,I])cccc1", "[CH1](=O)", "[#8][#8]",
  "[C;!R]=[N;!R]", "[N!R]=[N!R]", "[#6](=O)[#6](=O)", "[#16][#16]", "[#7][NH2]", "C(=O)N[NH2]",
  "[#6]=S", "[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]",
  "C1(=[O,N])C=CC(=[O,N])C=C1", "C1(=[O,N])C(=[O,N])C=CC=C1", "a21aa3a(aa1aaaa2)aaaa3",
  "a31a(a2a(aa1)aaaa2)aaaa3", "a1aa2a3a(a1)A=AA=A3=AA=A2", "c1cc([NH2])ccc1",
  "[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]",
  "I", "OS(=O)(=O)[O-]", "[N+](=O)[O-]", "C(=O)N[OH]",
  "C1NC(=O)NC(=O)1", "[SH]", "[S-]", "c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]",
  "c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]", "[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1",
  "[CR1]1[CR1][CR1]cc[CR1][CR1]1", "[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1",
  "[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1", "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
  "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1", "C#C",
  "[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]", "[$([N+R]),$([n+R]),$([N+]=C)][O-]",
  "[#6]=N[OH]", "[#6]=NOC=O", "[#6](=O)[CX4,CR0X3,O][#6](=O)", "c1ccc2c(c1)ccc(=O)o2",
  "[O+,o+,S+,s+]", "N=C=O", "[NX3,NX4][F,Cl,Br,I]", "c1ccccc1OC(=O)[#6]", "[CR0]=[CR0][CR0]=[CR0]",
  "[C+,c+,C-,c-]", "N=[N+]=[N-]", "C12C(NC(N1)=O)CSC2", "c1c([OH])c([OH,NH2,NH])ccc1", "P",
  "[N,O,S]C#N", "C=C=O", "[Si][F,Cl,Br,I]", "[SX2]O", "[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)",
  "O1CCCCC1OC2CCC3CCCCC3C2", "N=[CR0][N,n,O,S]",
  "[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2",
  "C=[C!r]C#N", "[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1",
  "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1",
  "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])", "[OH]c1ccc([OH,NH2,NH])cc1",
  "c1ccccc1OC(=O)O", "[SX2H0][N]", "c12ccccc1(SC(S)=N2)", "c12ccccc1(SC(=S)N2)", "c1nnnn1C=O",
  "s1c(S)nnc1NC=O", "S1C=CSC1=S", "C(=O)Onnn", "OS(=O)(=O)C(F)(F)F", "N#CC[OH]", "N#CC(=O)",
  "S(=O)(=O)C#N", "N[CH2]C#N", "C1(=O)NCC1", "S(=O)(=O)[O-,OH]", "NC[F,Cl,Br,I]", "C=[C!r]O",
  "[NX2+0]=[O+0]", "[OR0,NR0][OR0,NR0]", "C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]", "[CX2R0][NX3R0]",
  "c1ccccc1[C;!R]=[C;!R]c2ccccc2", "[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]",
  "[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]",
  "[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]", "[*]=[N+]=[*]", "[SX3](=O)[O-,OH]",
  "N#N", "F.F.F.F", "[R0;D2][R0;D2][R0;D2][R0;D2]", "[cR,CR]~C(=O)NC(=O)~[cR,CR]", "C=!@CC=[O,S]",
  "[#6,#8,#16][#6](=O)O[#6]", "c[C;R0](=[O,S])[#6]", "c[SX2][C;!R]", "C=C=C",
  "c1nc([F,Cl,Br,I,S])ncc1", "c1ncnc([F,Cl,Br,I,S])c1", "c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])",
  "[#6]S(=O)(=O)c1ccc(cc1)F", "[15N]", "[13C]", "[18O]", "[34S]"
  };

  // clang-format on

  QEDproperties WEIGHT_MAX =
      QEDproperties(0.50, 0.25, 0.00, 0.50, 0.00, 0.50, 0.25, 1.00);
  // currently only implemented with weight_mean
  QEDproperties WEIGHT_MEAN =
      QEDproperties(0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95);
  QEDproperties WEIGHT_NONE =
      QEDproperties(1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00);
  // Compute the RDKit molecule object for aliphatic rings
  std::unique_ptr<RDKit::RWMol> getAliphaticRingsMol();
  // Note for getStructurealAlertMols and getAcceptorMols:
  // Store the RDKit molecule itself and not the pointer because using a pointer
  // requires a pointer dereference on each iteration of the loop in
  // calcProperties whereas storing the molecule itself in the vector should
  // have the molecules contiguously in memory and better cache locality.
  // Also if using a unique_ptr, code working with the unique_ptr would need to
  // handle ownership and make sure it is moved correctly in the loop.
  //
  // Compute the RDKit molecule objects for the structural alerts
  std::vector<std::unique_ptr<RDKit::RWMol>> getStructuralAlertMols();
  // Compute the RDKit molecule objects for hydrogen bond acceptors
  std::vector<std::unique_ptr<RDKit::RWMol>> getAcceptorMols();
  // Calculate the properties needed for the QED descriptor
  QEDproperties calcProperties(const RDKit::ROMol &mol,
                               std::vector<RDKit::RWMol> acceptorMols,
                               std::unique_ptr<RDKit::RWMol> aliphaticRingMol,
                               std::vector<RDKit::RWMol> alertMols);
  // Compute the asymmetric double sigmoidal function using the value of the
  // descriptor of interest (the parameter `x` in the function) and the
  // adsParameters for that descriptor of interest
  double calcADS(float x, std::string adsParameterKey);
};
} // namespace duckdb_rdkit
