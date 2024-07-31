.open hundredkmol.duckdb 
.timer on
.mode csv
.output normal_substruct_mol.csv
LOAD duckdb_rdkit;
SET enable_progress_bar = false;


SELECT is_substruct(rdkit_mol, 'O=CNCCc1ccccc1') as s, count(*) from molecule group by s;
SELECT is_substruct(rdkit_mol, 'C1CCCC1') as s, count(*) from molecule group by s;
SELECT is_substruct(rdkit_mol, 'CC(C)C[C@H](NC(=O)[C@@H]1CCCN1C(=O)[C@@H]([NH3+])C(C)C)C([O-])=O') as s, count(*) from molecule group by s;
SELECT is_substruct(rdkit_mol, 'O=C1OC2=C(C=CC=C2)C=C1') as s, count(*) from molecule group by s;
SELECT is_substruct(rdkit_mol, 'OC(=O)C1=CC=CC=C1O') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1CCC2CCCC2C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1OC2=C(O1)C=CC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CCNC') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'ClC1=CC=CC=C1Cl') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'COc1cc2c(Nc3cc(CC(=O)Nc4cccc(F)c4F)[nH]n3)ncnc2cc1OCCCN(CCO)CC(C)C') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'FC(F)(F)C1=CC=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N1C=CC2=C1N=CC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N1N=CC2=CC=CC=C12') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'NC1=C2N=CNC2=NC=N1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'NC1=CC=NC(N)=N1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O=C1NC2=NC=NC=C2C=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'OCCCC1=CNC2=CC=CC=C12') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O=CNC1=NC=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N1C=CC=N1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'ONC') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1CCCC1.C2=CC=CC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1=CN(C=N1)C2=CC=CC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CBr') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CN1CCN(CC1)Cc2ccc(cc2)C(=O)Nc3ccc(C)c(Nc4nccc(n4)-c5cccnc5)c3') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CN(C(C)=O)C1=CC=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C[S]') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'NC1CCCCC1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'NC1=CC=CN=C1N') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'OC1CCCNC1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O=C1NC2=CC=CC=C2\C1=C\C3=CC=CN3') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O=C1OCCN1C2=CC=CC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'OC(=O)C1=CC=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CCNCCC1=CC=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O=C1NCCCN1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1=CC2=C(C=C1)C=NC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1=CC2=C(C=C1)N=CC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N=C=S') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1CC2=CC=CC=C2O1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1COC2=CC=CC=C2C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CC12CCC3C(CCC4=C3C=CC(O)=C4)C1CCC2O') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CCCC1=CC=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CCl') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CC(=O)NC1=CC=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N1N=CC2=C1C=CC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, '*NC') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'SC1=NC=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CNCC') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O=C1NCCCCN1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1C2=CC=CC=C2C3=C1C=NN3') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1CC1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, '[H]N') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O1C=CC2=CC=CC=C12') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'OCC') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1CCCCC1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1CC2=CC=CC=C2C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1=CN=CN=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CN1C=CC=N1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O=C1NC2=CC=CN=C2N1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N1C=NC2=C1C=CC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O=[N+](O)c1ccccc1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CCCC') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N1C=CC=C1.C2=CC=CC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N1C=NC2=CC=CC=C12') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O=C') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'OC(=O)C1=C(O)C=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'NC1=NC=CC=N1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O=C1NC2=CC=CC=C2N1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N1C=CC2=CC=CC=C12') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'S1C=NC2=CC=CC=C12') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CCC') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'FC') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CN1C=CC2=C1C=CC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1CNCCN1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1CCC2=CC=CC=C2C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1CCNC1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N1C=CC2=C1C=CC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'BrC1=CC=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CI') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'OC') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'S1C2=CC=CC=C2N=C1C1=CC=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CN') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, '[H]C(=C([H])C1=CC(O)=CC(O)=C1)C2=CC=C(O)C=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O=S') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N1C=CC=C1.C2=CC3=C(C=C2)C=CC=C3') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, '[H]CN') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1NCC2=C1C=CC3=C2C4=CC=CC=C4N3') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'N1C=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'OC1=CC=C(C=CC2=CC(O)=CC(O)=C2)C=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1=CC=NC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CC12CCC3C(CCC4=CC(O)=CC=C34)C1CCC2O') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CNC') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'NC') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1=CC=CC=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'C1=CC2=C(C=C1)C=CC=C2') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'O=C1OC2=CC=CC=C2C=C1') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CC') as s, count(*) from molecule group by s;
-- SELECT is_substruct(rdkit_mol, 'CC(C)C1(C)SC(Nc2ccccc2C(F)(F)F)=NC1=O') as s, count(*) from molecule group by s;




-- chembl_33
SELECT count(*) FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN compound_properties cp ON cp.molregno=m.molregno
      WHERE is_substruct(m.rdkit_mol, 'O=CNCCc1ccccc1');




SELECT a.standard_type, avg(a.value), count(a.value), a.relation, m.rdkit_mol FROM molecule m
INNER JOIN activities a ON a.molregno=m.molregno
INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
INNER JOIN compound_properties cp ON cp.molregno=m.molregno
WHERE is_substruct(m.rdkit_mol,'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
GROUP BY m.rdkit_mol, a.relation, a.standard_type;

