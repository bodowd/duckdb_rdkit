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






