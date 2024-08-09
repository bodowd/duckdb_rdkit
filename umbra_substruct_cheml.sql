.timer on
.open um.duckdb


-- Query 1

-- umbra
SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.umbra_mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE umbra_is_substruct(m.umbra_mol, 'COc1cc2c(Nc3cc(CC(=O)Nc4cccc(F)c4F)[nH]n3)ncnc2cc1OCCCN(CCO)CC(C)C');

-- standard
SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.rdkit_mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE is_substruct(m.rdkit_mol, 'COc1cc2c(Nc3cc(CC(=O)Nc4cccc(F)c4F)[nH]n3)ncnc2cc1OCCCN(CCO)CC(C)C');

-- postgres
SELECT pbd.prediction_method, a.value, a.units, a.type, a.relation, m.rdkit_mol FROM compound_structures m
  INNER JOIN activities a ON a.molregno=m.molregno
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  WHERE m.rdkit_mol@='COc1cc2c(Nc3cc(CC(=O)Nc4cccc(F)c4F)[nH]n3)ncnc2cc1OCCCN(CCO)CC(C)C';

-- Query 2

-- umbra
SELECT count(*) FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      WHERE umbra_is_substruct(m.umbra_mol, 'O=CNCCc1ccccc1');

-- standard 
SELECT count(*) FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      WHERE is_substruct(m.rdkit_mol, 'O=CNCCc1ccccc1');

-- postgres
SELECT count(*) FROM compound_structures m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      WHERE m.rdkit_mol@>'O=CNCCc1ccccc1';


-- Query 3

-- umbra
SELECT avg(a.value), stddev(a.value), a.units,a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.umbra_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE umbra_is_substruct(m.umbra_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
      GROUP BY m.umbra_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

-- normal
SELECT avg(a.value), stddev(a.value), a.units,a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE is_substruct(m.rdkit_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;


-- postgres
SELECT avg(a.value), stddev(a.value), a.units,a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM compound_structures m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE m.rdkit_mol@>'CC(=O)Nc1nnc(S(N)(=O)=O)s1'
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

-- Query 4

-- umbra
SELECT avg(a.value), stddev(a.value), a.units,a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.umbra_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE umbra_is_substruct(m.umbra_mol, 'O=CNCCc1ccccc1')
      GROUP BY m.umbra_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

-- normal
SELECT avg(a.value), stddev(a.value), a.units,a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE is_substruct(m.rdkit_mol, 'O=CNCCc1ccccc1')
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;


-- postgres
SELECT avg(a.value), stddev(a.value), a.units,a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM compound_structures m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE m.rdkit_mol@>'O=CNCCc1ccccc1'
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;


-- Query 5
-- umbra
  SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.umbra_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE umbra_is_substruct(m.umbra_mol, 'N1C=CC=N1')
      GROUP BY m.umbra_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

-- standard
  SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM molecule m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE is_substruct(m.rdkit_mol, 'N1C=CC=N1')
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;
  

-- postgres
SELECT avg(a.value), stddev(a.value), a.units, a.type, count(a.value), a.relation, bs.site_name, ys.assay_organism, m.rdkit_mol FROM compound_structures m
      INNER JOIN activities a ON a.molregno=m.molregno
      INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
      INNER JOIN binding_sites bs ON pbd.site_id=bs.tid
      INNER JOIN assays ys ON ys.tid=bs.tid
      WHERE m.rdkit_mol@>'N1C=CC=N1'
      GROUP BY m.rdkit_mol, a.relation, a.units, a.type, bs.site_name, ys.assay_organism;

