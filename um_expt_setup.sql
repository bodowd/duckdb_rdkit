.open um.duckdb
.timer on;
install postgres;
load postgres;

ATTACH 'dbname=chembl_33 host=localhost user=postgres password=example port=5435' AS pg (TYPE POSTGRES, READ_ONLY); 

COPY (SELECT * FROM pg.compound_structures) TO 'compound_structures.parquet' (FORMAT parquet);
COPY (SELECT * FROM pg.binding_sites) TO 'binding_sites.parquet' (FORMAT parquet);
COPY (SELECT * FROM pg.activities) TO 'activities.parquet' (FORMAT parquet);
COPY (SELECT * FROM pg.assays) TO 'assays.parquet' (FORMAT parquet);
COPY (SELECT * FROM pg.predicted_binding_domains) TO 'predicted_binding_domains.parquet' (FORMAT parquet);


create table molecule as select * from 'compound_structures.parquet';
create table binding_sites as select * from 'binding_sites.parquet';
create table activities as select * from 'activities.parquet';
create table assays as select * from 'assays.parquet';
create table predicted_binding_domains as select * from 'predicted_binding_domains.parquet';


alter table molecule add column umbra_mol UmbraMol;
alter table molecule drop column rdkit_mol;
alter table molecule add column rdkit_mol Mol;

update molecule set umbra_mol=umbra_mol_from_smiles(canonical_smiles);
update molecule set rdkit_mol=mol_from_smiles(canonical_smiles);
