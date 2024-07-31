.timer on
.open chembl_33.duckdb

--Query1
select molregno,canonical_smiles, rdkit_mol, umbra_mol from molecule where umbra_is_exact_match(umbra_mol,'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');
--Query2

select molregno,canonical_smiles, rdkit_mol, umbra_mol from molecule where umbra_is_exact_match(umbra_mol,'CCC');
--Query3

SELECT pbd.prediction_method, a.value, a.relation, m.umbra_mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  INNER JOIN compound_properties cp ON cp.molregno=m.molregno
  WHERE umbra_is_exact_match(m.umbra_mol, 'COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC');
--Query4

SELECT avg(a.value), count(a.value), a.relation, m.umbra_mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  INNER JOIN compound_properties cp ON cp.molregno=m.molregno
  WHERE umbra_is_exact_match(m.umbra_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
  GROUP BY m.umbra_mol, a.relation;






-- normal exact match
select molregno, canonical_smiles, rdkit_mol, umbra_mol from molecule where is_exact_match(mol,'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');

select molregno, canonical_smiles, rdkit_mol, umbra_mol from molecule where is_exact_match(mol,'CCC');




SELECT pbd.prediction_method, a.value, a.relation, m.rdkit_mol FROM molecule m
    INNER JOIN activities a ON a.molregno=m.molregno
    INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
    INNER JOIN compound_properties cp ON cp.molregno=m.molregno
    WHERE is_exact_match(m.rdkit_mol, 'COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC');

SELECT avg(a.value), count(a.value), a.relation, m.rdkit_mol FROM molecule m
INNER JOIN activities a ON a.molregno=m.molregno
INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
INNER JOIN compound_properties cp ON cp.molregno=m.molregno
WHERE is_exact_match(m.rdkit_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
GROUP BY m.rdkit_mol, a.relation;



-- postgres
SELECT a.standard_type, avg(a.value), count(a.value), a.relation, m.rdkit_mol FROM compound_structures m
INNER JOIN activities a ON a.molregno=m.molregno
INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
INNER JOIN compound_properties cp ON cp.molregno=m.molregno
WHERE m.rdkit_mol@>'CC(=O)Nc1nnc(S(N)(=O)=O)s1'
GROUP BY m.rdkit_mol, a.relation, a.standard_type;


