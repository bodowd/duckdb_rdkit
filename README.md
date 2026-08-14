# An RDKit extension for DuckDB

[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)

## Getting started

- Download [duckdb](https://duckdb.org/install/?platform=linux&environment=cli)
- Start the duckdb CLI in your terminal and run: 
```shell
> install duckdb_rdkit from community;
> load duckdb_rdkit;
> select is_substruct('CC', 'C');
┌─────────────────────────┐
│ is_substruct('CC', 'C') │
│         boolean         │
├─────────────────────────┤
│ true                    │
└─────────────────────────┘
```

> [!NOTE]
> Currently supported platforms are osx_arm64, linux_arm64, and linux_amd64

---

This extension integrates RDKit into DuckDB to enable you to do
cheminformatics work with DuckDB.

## Currently supported functionality:

### Searches

- `is_exact_match(mol1, mol2)`: exact structure search. Returns true if the two molecules are the same. (Chirality sensitive search is not on)
  - Note: if you are looking for very specific capabilities with exact match with regards
    to stereochemistry or tautomers, the `RegistrationHash` (https://rdkit.org/docs/source/rdkit.Chem.RegistrationHash.html)
    might be an option to consider. You would need to write this to your DB and
    then you can do a simple VARCHAR based search on those columns.
- `is_substruct(mol1, mol2)`: returns true if mol2 is a substructure of mol1.

### File formats

#### SDF

- There are two ways to query `.sdf` files with SQL.
  These can be used to extract, transform, and load data
  into a duckdb file for faster subsequent queries, or to directly query the
  sdf to explore the data.
  - `read_sdf(path/to/sdf/file, COLUMNS={column_name: LogicalType});`
    Using the `read_sdf` function, the properties of interest in the sdf file
    can be explicitly defined. If a record does not have the specified property,
    a null value will be returned. The `'Mol'` type will indicate to the
    extension that the molecules in the records should be extracted and returned.
    - Example: `SELECT * FROM read_sdf(path/to/file, COLUMNS={desired_col: 'VARCHAR', mol: 'Mol'});`

  - Automatic detection of `sdf` files. This will execute the query against
    the sdf file when the extension `.sdf` is detected.

    In this case, the extension
    will guess what the schema is. If the schema is not homogeneous, it is possible
    that the automatic detection will miss certain properties in the SDF.

    In this case, it is better to use the `read_sdf` function in order to make
    sure the property of interest is extracted. This is not a problem if the
    schema is uniform throughout the sdf file.

    The molecule column is named `mol`.

  - Example: `SELECT mol, id FROM 'test.sdf';`


### Types

- `Mol`: the internal duckdb_rdkit representation of a RDKit molecule.
  - Currently only SMILES can be converted to `Mol`. This can be done with
    `mol_from_smiles`, or by casts (i.e. inserting a SMILES string into a
    column that expects `Mol` or `'CC::mol'`).

> [!IMPORTANT]  
> The duckdb_rdkit molecule representation has additional metadata and cannot
> be read directly by RDKit. You will get an error. You can use `mol_to_rdkit_mol`
> to convert the duckdb_rdkit molecule representation into one that is RDKit compatible.

### Molecule conversion functions

- `mol_from_smiles(SMILES)`: returns a molecule for a SMILES string. Returns NULL if mol cannot be made from SMILES
- `mol_to_smiles(mol)`: returns the SMILES string for a RDKit molecule
- `mol_to_rdkit_mol(mol)`: returns the binary RDKit molecule in hexadecimal representation
  - duckdb_rdkit has its own binary representation of molecules, which differs from RDKit’s format.
    Use this function to extract a molecule from duckdb_rdkit and convert it
    into a format compatible with RDKit. The returned value can be passed
    to RDKit's `Chem.Mol` function for further processing in Python.

### Molecule descriptors

- `mol_logp(mol)`: returns the Wildman-Crippen LogP estimate for a molecule
- `mol_exactmw(mol)`: returns the exact molecular weight
- `mol_amw(mol)`: returns the approximate molecular weight
- `mol_tpsa(mol)`: returns the topological polar surface area
- `mol_hba(mol)`: returns the number of H-bonds acceptors
- `mol_hbd(mol)`: returns the number of H-bonds donors
- `mol_num_rotatable_bonds(mol)`: returns the number of rotatable bonds
- `mol_qed(mol)`: returns the quantitative estimate of drug-likeness (QED) of the molecule
  - currently only implements the "mean weight" of the ADS parameters from the paper Quantifying the chemical beauty of drugs by Bickerton, et al.


### Building duckdb_rdkit

First, clone this repository with recurse submodules to pull duckdb and the
extension-ci-tools repositories

```shell
> git clone --recurse-submodules https://github.com/bodowd/duckdb_rdkit.git
```

Then run:

```shell
> GEN=ninja make release
```
This will build RDKit and statically link it to duckdb.

For further information on building duckdb from source,
you can visit https://duckdb.org/docs/dev/building/overview.html

### In the python client

Then test it out:

```python
import duckdb

con = duckdb.connect()
con.execute("INSTALL duckdb_rdkit FROM community;")
con.execute("LOAD duckdb_rdkit;")
# should return true
print(con.sql("SELECT is_exact_match('C', 'C');"))
# should return false
print(con.sql("SELECT is_exact_match('C', 'CO');"))

```

## Running the tests

Different tests can be created for DuckDB extensions. The primary way of testing DuckDB extensions should be the SQL tests in `./test/sql`. These SQL tests can be run using:

```sh
make test
```
