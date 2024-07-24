.open demo.duckdb;
.timer on;
install postgres;
load postgres;
create table molecule (id integer primary key, smiles varchar);
insert into molecule select molregno, canonical_smiles from postgres_scan('dbname=chembl_33 host=localhost user=postgres password=example port=5435', 'public', 'compound_structures');

-- create the rdkit mol column 
alter table molecule add column mol Mol;
update molecule set mol=mol_from_smiles(smiles);

-- create the umbra mol column
ALTER TABLE molecule ADD COLUMN umbra_mol UmbraMol;
update molecule set umbra_mol=umbra_mol_from_smiles(mol);


-- search using rdkit mol object
select * from molecule where is_exact_match(mol,'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');
select * from molecule where is_exact_match(mol,'CCC');

-- using umbra mol 
select * from molecule where umbra_is_exact_match(umbra_mol,'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O');
select * from molecule where umbra_is_exact_match(umbra_mol,'CCC');




-- copy data from pg into duckdb database using these statements, but just change the table names
copy (select * from pg.binding_sites) to 'binding_sites.parquet' (format parquet);
create table binding_sites as select * from 'binding_sites.parquet';





-- joining three tables with an exact match search

SELECT pbd.prediction_method, a.value, a.relation, m.mol FROM molecule m
    INNER JOIN activities a ON a.molregno=m.id
    INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
    INNER JOIN compound_properties cp ON cp.molregno=m.id
    WHERE is_exact_match(m.mol, 'COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC');

SELECT pbd.prediction_method, a.value, a.relation, m.umbra_mol FROM molecule m
    INNER JOIN activities a ON a.molregno=m.id
    INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
    INNER JOIN compound_properties cp ON cp.molregno=m.id
    WHERE umbra_is_exact_match(m.umbra_mol, 'COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC');

-- same query but in postgres
SELECT pbd.prediction_method, a.value, a.relation, m.rdkit_mol FROM compound_structures m
    INNER JOIN activities a ON a.molregno=m.molregno
    INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
    INNER JOIN compound_properties cp ON cp.molregno=m.molregno
    WHERE m.rdkit_mol@='COc1cc(/C=C/C(=O)OCCCCCCN(C)CCCCOC(=O)c2c3ccccc3cc3ccccc23)cc(OC)c1OC';



SELECT pbd.prediction_method, a.value, a.relation, m.mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.id
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  INNER JOIN compound_properties cp ON cp.molregno=m.id
  WHERE is_exact_match(m.mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1');

-- with umbra mol
SELECT pbd.prediction_method, a.value, a.relation, m.umbra_mol FROM molecule m
  INNER JOIN activities a ON a.molregno=m.id
  INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
  INNER JOIN compound_properties cp ON cp.molregno=m.id
  WHERE umbra_is_exact_match(m.umbra_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1');




-- standard mol
SELECT avg(a.value), count(a.value), a.relation, m.mol FROM molecule m
INNER JOIN activities a ON a.molregno=m.id
INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
INNER JOIN compound_properties cp ON cp.molregno=m.id
WHERE is_exact_match(m.mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
GROUP BY m.mol, a.relation;

-- now with umbra mol
SELECT avg(a.value), count(a.value), a.relation, m.umbra_mol FROM molecule m
INNER JOIN activities a ON a.molregno=m.id
INNER JOIN predicted_binding_domains pbd ON pbd.activity_id=a.activity_id
INNER JOIN compound_properties cp ON cp.molregno=m.id
WHERE umbra_is_exact_match(m.umbra_mol, 'CC(=O)Nc1nnc(S(N)(=O)=O)s1')
GROUP BY m.umbra_mol, a.relation;
