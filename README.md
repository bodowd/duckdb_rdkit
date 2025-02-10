# An RDKit extension for DuckDB

[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)

This repository is based on https://github.com/duckdb/extension-template, check it out if you want to build and ship your own DuckDB extension.

---

This extension, duckdb_rdkit, integrates RDKit into DuckDB to enable you to do
cheminformatics work with DuckDB.

## Currently supported functionality:

### Types

- `Mol`: the internal duckdb_rdkit representation of a RDKit molecule.
  - Currently only SMILES can be converted to `Mol`. This can be done with
    `mol_from_smiles`, or by casts (i.e. inserting a SMILES string into a
    column that expects `Mol` or `'CC::mol'`).

> [!IMPORTANT]  
> The duckdb_rdkit molecule representation has additional metadata and cannot
> be read directly by RDKit. You will get an error. You can use `mol_to_rdkit_mol`
> to convert the duckdb_rdkit molecule representation into one that is RDKit compatible.

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

### Searches

- `is_exact_match(mol1, mol2)`: exact structure search. Returns true if the two molecules are the same. (Chirality sensitive search is not on)
  - Note: if you are looking for very specific capabilities with exact match with regards
    to stereochemistry or tautomers, the `RegistrationHash` (https://rdkit.org/docs/source/rdkit.Chem.RegistrationHash.html)
    might be an option to consider. You would need to write this to your DB and
    then you can do a simple VARCHAR based search on those columns.
- `is_substruct(mol1, mol2)`: returns true if mol2 is a substructure of mol1.

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
- `mol_qed(mol)`: returns the quantitative estimate of drug-likeness (QED) of the molecule
  - currently only implements the "mean weight" of the ADS parameters from the paper Quantifying the chemical beauty of drugs by Bickerton, et al.

## Getting started

Unfortunately, I haven't been able to find a way to make installing the duckdb_rdkit
extension as easy as `INSTALL` and `LOAD` as other duckdb
extensions may be.

I have only been able to successfully compile and test the extension on `linux_amd64`
and `osx_arm64`.

You can download the binaries from the releases, or build the extension from source.
The compiled binary built for the duckdb_rdkit extension is not signed. You may get
a warning about running unverified applications from the OS.

- [Building](#building) section.

- [Running](#running) the extension section.

## <a name="building"></a>Building

### Building duckdb_rdkit

First, clone this repository with recurse submodules to pull duckdb and the
extension-ci-tools repositories

```shell
git clone --recurse-submodules https://github.com/bodowd/duckdb_rdkit.git
```

To build the extension, you need to have RDKit installed.
The instructions below are derived from this post on the RDKit [blog](https://greglandrum.github.io/rdkit-blog/posts/2021-07-24-setting-up-a-cxx-dev-env.html).
The easiest way to install RDKit is with conda, and I used [miniforge](https://github.com/conda-forge/miniforge).

After installing conda, you can create a new
conda environment and then install the packages needed.
`linux_conda_env.yml` or `osx_conda_env.yml` can be used to create a conda
environment for building the extension.

```shell
# activate your conda env and then in your conda env run:
conda create -n rdkit_dev
conda activate rdkit_dev
# or use the osx_conda_env.yml if you are on osx
conda env update -f linux_conda_env.yml
```

After installing the prerequisite software, you can run:

```shell
GEN=ninja make
```

This will compile duckdb and the extension and you will find it in
the `build` folder.

For further information on building duckdb from source,
you can visit https://duckdb.org/docs/dev/building/overview.html

## <a name="running"></a> Running the extension

### In the CLI

If you want to run the duckdb binary you built from source from this
duckdb_rdkit extension repository, you can just run `./build/release/duckdb`.
This will already have the extension loaded in.

If you downloaded the compiled binaries from here, you will need to tell
duckdb where to find the RDKit shared object files. Otherwise, you may see errors like this:
`./duckdb: error while loading shared libraries: libRDKitDescriptors.so.1: cannot open shared object file: No such file or directory`

If you have your conda env activated:

```shell
# LINUX
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
# OSX
export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib:$DYLD_LIBRARY_PATH
```

If you don't have your conda env activated, you will need to find where
your installation has placed these files. For example, in `~/miniforge3/envs/my_rdkit_env/lib`.
You will need to add your path to `LD_LIBRARY_PATH` on Linux, or `DYLD_LIBRARY_PATH` on osx.

If you want to run with a different binary that does not have the extension already
installed and loaded, but rather point to this extension,
you'll need to tell duckdb where the extension is, and you also need to tell
it to run unsigned extensions.

Warning: I was not able to get the extension to run on the linux CLI binary downloaded
from duckdb's website. That seems to have been compiled for `linux_amd64_gcc4`,
and I was not successful compiling the extension for that.

Run duckdb with the `unsigned` flag on to run unsigned extensions.
More information here: https://duckdb.org/docs/extensions/overview.html#unsigned-extensions

```shell
duckdb -unsigned
```

Then load the extension with the path to the duckdb_extension file:

```shell
LOAD 'path/to/duckdb_rdkit.duckdb_extension'
```

Now confirm if the extension is working:

```shell
# should return true
SELECT is_exact_match('C', 'C');

# should return false
SELECT is_exact_match('C', 'CO');
```

### In the python client

Warning: On Linux, I was unable to get the client I installed via pip to load the
extension because it only seems to support loading extensions compiled for `linux_amd64_gcc4`.
I was able to get it loaded in duckdb installed via conda though. See [duckdb's website](https://duckdb.org/docs/api/python/overview.html#:~:text=The%20DuckDB%20Python%20API%20can,requires%20Python%203.7%20or%20newer.) for
more information.

See the duckdb [documentation](https://duckdb.org/docs/api/python/overview.html#:~:text=The%20DuckDB%20Python%20API%20can,requires%20Python%203.7%20or%20newer.)
for instructions on installing the python client.

You may need to tell duckdb where to find the RDKit shared object files.

```shell
# LINUX
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
# OSX
export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib:$DYLD_LIBRARY_PATH
```

Then test it out:

```python
import duckdb
con = duckdb.connect(config = {"allow_unsigned_extensions": "true"})
con.install_extension('/path/to/duckdb_rdkit.duckdb_extension')
con.load_extension('/path/to/duckdb_rdkit.duckdb_extension')
# should return true
con.sql("SELECT is_exact_match('C', 'C');")
# should return false
con.sql("SELECT is_exact_match('C', 'CO');")

```

## Running the tests

Different tests can be created for DuckDB extensions. The primary way of testing DuckDB extensions should be the SQL tests in `./test/sql`. These SQL tests can be run using:

```sh
make test
```
