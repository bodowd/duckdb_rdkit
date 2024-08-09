.timer on
.open um.duckdb

-- Query 1
-- umbra
select molregno,canonical_smiles, rdkit_mol, umbra_mol from molecule where umbra_is_exact_match(umbra_mol,'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');

-- standard
select molregno,canonical_smiles, rdkit_mol, umbra_mol from molecule where is_exact_match(rdkit_mol,'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');

-- postgres
select molregno,canonical_smiles, rdkit_mol from compound_structures where rdkit_mol@='Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');



-- Query 2
select molregno,canonical_smiles, rdkit_mol, umbra_mol from molecule where umbra_is_exact_match(umbra_mol,'CCC');


-- Query 3
-- umbra
SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.umbra_mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE umbra_is_exact_match(m.umbra_mol, 'COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC');

-- standard
SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.rdkit_mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE is_exact_match(m.rdkit_mol, 'COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC');

-- postgres
SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.rdkit_mol FROM compound_structures m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE m.rdkit_mol@='COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC';

-- Query 4
-- umbra
  SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.umbra_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE umbra_is_exact_match(m.umbra_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
      GROUP BY m.umbra_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

-- standard
SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE is_exact_match(m.rdkit_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;


-- postgres
SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM compound_structures m
  INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
  WHERE m.rdkit_mol@='CC(=O)Nc1nnc(S(N)(=O)=O)s1'
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;
