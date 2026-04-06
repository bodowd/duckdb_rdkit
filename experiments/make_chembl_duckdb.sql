INSTALL sqlite;
LOAD sqlite;

ATTACH '/home/bing/.data/chembl/35/chembl_35.db' as sqlite (TYPE sqlite, READ_ONLY);

CREATE TABLE compound_structures AS 
  SELECT molregno, 
    canonical_smiles, 
    mol_from_smiles(canonical_smiles) AS mol
  FROM sqlite.compound_structures;


CREATE TABLE activities AS
SELECT * FROM sqlite.activities;

CREATE TABLE predicted_binding_domains AS
SELECT * FROM sqlite.predicted_binding_domains;

CREATE TABLE binding_sites AS
SELECT * FROM sqlite.binding_sites;

