# An RDKit extension for DuckDB

[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)

This repository is based on https://github.com/duckdb/extension-template, check it out if you want to build and ship your own DuckDB extension.

---

This extension, duckdb_rdkit, allows you to use RDKit functionality within DuckDB.

## Currently supported functionality:

### Types

- Mol: an RDKit molecule. Currently, can only be created from a SMILES in a variety of ways: inserting a valid SMILES
  string into a column that expects Mol, type conversion such as 'CC'::mol, or the mol_from_smiles function.

### Searches

- is_exact_match(mol1, mol2) : exact structure search. Returns true if the two molecules are the same. (Chirality sensitive search is not on)
  - Note: if you are looking for very specific capabilities with exact match with regards
    to stereochemistry or tautomers, the `RegistrationHash` (https://rdkit.org/docs/source/rdkit.Chem.RegistrationHash.html)
    might be an option to consider. You would need to write this to your DB and
    then you can do a simple VARCHAR based search on those columns.
- is_substruct(mol1, mol2): returns true if mol2 is a substructure of mol1.

### Molecule conversion functions

- mol_from_smiles(SMILES): returns a molecule for a SMILES string. Returns NULL if mol cannot be made from SMILES
- mol_to_smiles(mol): returns the SMILES string for a RDKit molecule

## Some important information on getting started

Unfortunately, I haven't been able to find a way to make installing the duckdb_rdkit
extension as easy as `INSTALL` and `LOAD` as other duckdb
extensions may be.

Long story short, I suggest building the extension from source,
and you can find more information in the
[Building](#building) section.

### Linux

The duckdb binary from duckdb's website, seems to only support
loading extensions that were compiled with gcc4 (`linux_amd64_gcc4` platform in my
case). I was only able to get the extension
to compile with newer versions of GCC, possibly due to the RDKit dependency.
Therfore, I was not able to load the duckdb_rdkit extension in a duckdb binary that
I downloaded from their website.

Probably the best way is to build the extension from source. This will give you
a duckdb binary, as well as the duckdb_rdkit.duckdb_extension binary. The compilation
will consider your OS.

I will also work on making the compiled binaries from the github actions
workflows available for certain platforms.

### MacOS

I suggest building from source.

The compiled binary for the duckdb rdkit extension is not signed. You may get
a warning about running unverified applications from the OS.

Currently, I am unable to test if downloading the duckdb binary from their website can
load the duckdb_rdkit extension.

I will also work on making the compiled binaries from the github actions
workflows available for certain osx platforms.

### Windows

I have not tested anything with Windows yet.

## <a name="building"></a>Building

### Building duckdb_rdkit

First, clone this repository with recurse submodules to pull duckdb and the
extension-ci-tools repositories

```shell
git clone --recurse-submodules https://github.com/bodowd/duckdb_rdkit.git
```

To build the extension, you need to have RDKit installed.
The instructions below are derived from this post on the RDKit [blog].
The easiest way to install RDKit is with conda.
I used miniforge to install things: https://github.com/conda-forge/miniforge

After installing conda, you can create a new
conda environment and then install the packages needed.
As of August 2024, I found installing these packages worked (librdkit-dev seems to have the relevant header files).

```shell
# activate your conda env and then in your conda env run:
conda install -c conda-forge -y boost-cpp boost cmake rdkit eigen librdkit-dev
```

After installing the prerequisite software, you can run:

```shell

GEN=ninja make
```

This will compile duckdb and the extension and you will find it in
the `build` folder.

For further information on building duckdb from source,
you can visit https://duckdb.org/docs/dev/building/overview.html

## Running the extension

### In the CLI

If you want to run the duckdb binary you built from source from this
duckdb_rdkit extension repository, you can just run `./build/release/duckdb`.
This will already have the extension loaded in.

If you downloaded the compiled binaries from here, you will need to tell
duckdb where to find the RDKit shared object files. Otherwise, you may see this error:
`./duckdb: error while loading shared libraries: libRDKitDescriptors.so.1: cannot open shared object file: No such file or directory`

If you have your conda env activated:

```shell
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
```

If you don't have your conda env activated, you will need to find where
your installation has placed these files. For example, in `~/miniforge3/envs/my_rdkit_env/lib`

If you want to run with a different binary, but point to this extension,
you'll need to tell duckdb where the extension is, and you also need to tell
it to run unsigned extensions.

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
extension because it only seems to support loading extensions compiled with gcc 4.
I was able to get the extension to load in the python client on OSX though.

See the duckdb documentation for instructions on installing the python client.

You may need to tell duckdb where to find the RDKit shared object files.

```shell
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
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

<!-- #### Issue with building on MacOS 14.3 -->
<!---->
<!-- There was an issue building on MacOS 14.3 where a header from boost could not be found. -->
<!-- If you have trouble, you can try creating a conda env using the `starter_conda_env.yml` (from DavidACosgrove in the RDKit discussions). -->
<!-- This includes boost libraries that are needed by RDKit. Then install RDKit into that env. And then run `make` in -->
<!-- that env. -->
<!---->
<!-- [blog]: https://greglandrum.github.io/rdkit-blog/posts/2021-07-24-setting-up-a-cxx-dev-env.html -->

## Running the tests

Different tests can be created for DuckDB extensions. The primary way of testing DuckDB extensions should be the SQL tests in `./test/sql`. These SQL tests can be run using:

```sh
make test
```
