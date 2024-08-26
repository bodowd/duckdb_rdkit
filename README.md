# An RDKit extension for DuckDB

[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)

This repository is based on https://github.com/duckdb/extension-template, check it out if you want to build and ship your own DuckDB extension.

---

This extension, duckdb_rdkit, allows you to use RDKit functionality within DuckDB.

## Currently supported functionality:

### Types

- Mol: an rdkit molecule. Can be created from a SMILES in a variety of ways: inserting a valid SMILES
  string into a column that expects Mol, type conversion such as 'CC'::mol, or the mol_from_smiles function.

### Searches

- is_exact_match(mol1, mol2) : exact structure search. Returns true if the two molecules are the same. (Chirality sensitive search is not on)
- is_substruct(mol1, mol2)): returns true if mol2 is a substructure of mol1.

### Molecule conversion functions

- mol_from_smiles(SMILES): returns a molecule for a SMILES string. Returns NULL if mol cannot be made from SMILES
- mol_to_smiles(mol): returns the SMILES string for a rdkit molecule

## Building

### Building duckdb with RDKit extension

To build the extension, you need to have RDKit installed. The easiest way
to install RDKit is with conda. I used miniforge to install things: https://github.com/conda-forge/miniforge

After installing it and having conda working, you can create a new
conda environment and then install the packages needed. The instructions are derived from this post on the RDKit [blog].
As of August 2024, I am using this command to install the libraries needed:

```shell
# activate your conda env and then in your conda env run:
conda install -c conda-forge -y boost-cpp boost cmake rdkit eigen librdkit-dev
```

(librdkit-dev seems to have the relevant header files)

You can visit https://duckdb.org/docs/dev/building/overview.html
to find more information on how to build duckdb from source.

After installing the prerequisite software, you can run `GEN=ninja make`.
This will compile duckdb and the extension and you will find it in
the `build` folder.

### Running in the CLI

You can either run the duckdb executable you compiled found in
`build/release` and that will have the duckdb_rdkit extension already
installed and loaded, or you can run a different duckdb executable, and then
point that to the the duckdb_rdkit extension file.

If you want to run a different executable that you installed and then load in
the duckdb_rdkit extension, you can do the following:

Run duckdb with the `unsigned` flag on to run unsigned extensions.
More information here: https://duckdb.org/docs/extensions/overview.html#unsigned-extensions

```shell
duckdb -unsigned
```

Then load the extension with the path to the duckdb_extension file:

```shell
LOAD 'path/to/duckdb_rdkit.duckdb_extension'
```

Now try

```shell
# should return true
SELECT is_exact_match('C', 'C');

# should return false
SELECT is_exact_match('C', 'CO');
```

#### Running in python client:

See duckdb documentation for instructions on installing the python client.

You will need to tell duckdb where to find the RDKit shared object files.

<!-- ```shell -->
<!-- export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH -->
<!-- ``` -->

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

#### Running the extension by downloading the binary extension file, and not building it yourself

If you download the binary extension file, and do not build it yourself, you will need to
tell duckdb where the RDKit shared object libraries are.

When RDKit is installed with conda, they can be found in `$CONDA_PREFIX/lib`

```shell
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
```

Then you can run according to the instructions as described in the above sections on running
in the CLI or in the Python client.

#### Issue with building on MacOS 14.3

There was an issue building on MacOS 14.3 where a header from boost could not be found.
If you have trouble, you can try creating a conda env using the `starter_conda_env.yml` (from DavidACosgrove in the RDKit discussions).
This includes boost libraries that are needed by RDKit. Then install RDKit into that env. And then run `make` in
that env.

[blog]: https://greglandrum.github.io/rdkit-blog/posts/2021-07-24-setting-up-a-cxx-dev-env.html

## Running the tests

Different tests can be created for DuckDB extensions. The primary way of testing DuckDB extensions should be the SQL tests in `./test/sql`. These SQL tests can be run using:

```sh
make test
```
