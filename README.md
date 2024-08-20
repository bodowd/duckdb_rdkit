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

- is_exact_match(mol1, mol2) : exact structure search. Returns boolean

### Molecule conversion functions

- mol_from_smiles(SMILES): returns a molecule for a SMILES string. Returns NULL if mol cannot be made from SMILES
- mol_to_smiles(mol): returns the SMILES string for a rdkit molecule

## Building

### Building RDKit with static libraries on Linux

Create a conda environment:

```shell
# not sure if all of these are needed, like py-boost
conda create -n rdkit_dev -c conda-forge -y boost-cpp boost cmake cairo eigen libboost py-boost
```

Clone RDKit:

```shell
git clone https://github.com/rdkit/rdkit.git
```

Build and install RDKit:

```shell
cd rdkit
mkdir build
cd build
RUN cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DRDK_INSTALL_INTREE=ON \
  -DRDK_BUILD_PYTHON_WRAPPERS=OFF\
  -DRDK_BUILD_FUZZ_TARGETS=OFF \
  -DRDK_INSTALL_STATIC_LIBS=ON \
  -DBoost_USE_STATIC_LIBS=ON \
  -DRDK_BUILD_CPP_TESTS=OFF \
  -DBoost_NO_SYSTEM_PATHS=ON \
  -DCMAKE_INCLUDE_PATH="${CONDA_PREFIX}/include" \
  -DCMAKE_LIBRARY_PATH="${CONDA_PREFIX}/lib" \
  -DCMAKE_PREFIX_PATH=$CONDA_PREFIX

make install
```

After building RDKit, should be able to find static libraries:

```shell
find . -type f -name "*.a"
```

returns for example:

```shell
./Code/ChemicalFeatures/libRDKitChemicalFeatures_static.a
./Code/DataStructs/libRDKitDataStructs_static.a
./Code/GraphMol/MolHash/libRDKitMolHash_static.a
./Code/GraphMol/MolCatalog/libRDKitMolCatalog_static.a
./Code/GraphMol/FilterCatalog/libRDKitFilterCatalog_static.a
./Code/GraphMol/MolChemicalFeatures/libRDKitMolChemicalFeatures_static.a
...
```

### Building duckdb with RDKit extension dynamic linking

Tried building in the host in conda env (which has the RDKit shared objects
in ~/miniforge3/envs/rdkit_dev/lib) then running it in the docker conda, which has
the RDKit shared objects in /opt/conda/envs/rdkit_dev/lib.

```docker
FROM continuumio/miniconda3
RUN conda install -c conda-forge mamba && \
  conda create -n rdkit_dev -c conda-forge -y boost-cpp boost cmake rdkit=2023.09.4 eigen
```

Then build the container:

```shell
docker build -t duckdb_rdkit_image .
```

In the host, I activated the conda env there and built duckdb with extension.

Then run the docker container:

```shell
docker run -it --rm --name=duckdb_rdkit_image --mount type=bind,source=${PWD},target=/src duckdb_rdkit_image bash
```

In the container set LD_LIBRARY_PATH to the shared libraries. The container
has a different $CONDA_PREFIX

```shell
export LD_LIBRARY_PATH=$CONDA_PREFIX/envs/rdkit_dev/lib:$LD_LIBRARY_PATH
```

Then run duckdb:

```shell
cd src
./build/release/duckdb
```

The extension then should be working with `select mol_from_smiles('C')`

TODO: try again the other way around. Build in the container, then run
in the host

### Building duckdb with RDKit extension static linking

#### Install RDKit

The instructions are derived from this post on the RDKit [blog].

The easiest way is to install via conda/mamba.

Summary: Install conda and mamba. Create and activate a conda environment. Install RDkit in there.
Run `make` to build duckdb and the extension.

#### Issue with building on MacOS 14.3

There was an issue building on MacOS 14.3 where a header from boost could not be found.
You can try creating a conda env using the `starter_conda_env.yml` (from DavidACosgrove in the RDKit discussions).
This includes boost libraries that are needed by RDKit. Then install RDKit into that env. And then run `make` in
that env.

[blog]: https://greglandrum.github.io/rdkit-blog/posts/2021-07-24-setting-up-a-cxx-dev-env.html

#### Build duckdb with RDKit extension steps

To build the extension, run:

```sh
make
```

The main binaries that will be built are:

```sh
./build/release/duckdb
./build/release/test/unittest
./build/release/extension/duckdb_rdkit/duckdb_rdkit.duckdb_extension
```

- `duckdb` is the binary for the duckdb shell with the extension code automatically loaded.
- `unittest` is the test runner of duckdb. Again, the extension is already linked into the binary.
- `duckdb_rdkit.duckdb_extension` is the loadable binary as it would be distributed.

## Running the extension

To run the extension code, simply start the shell with `./build/release/duckdb`.

Now we can use the features from the extension directly in DuckDB.

```
D select mol_from_smiles('CC') as result;
```

## Running the tests

Different tests can be created for DuckDB extensions. The primary way of testing DuckDB extensions should be the SQL tests in `./test/sql`. These SQL tests can be run using:

```sh
make test
```

### Installing the deployed binaries

To install your extension binaries from S3, you will need to do two things. Firstly, DuckDB should be launched with the
`allow_unsigned_extensions` option set to true. How to set this will depend on the client you're using. Some examples:

CLI:

```shell
duckdb -unsigned
```

Python:

```python
con = duckdb.connect(':memory:', config={'allow_unsigned_extensions' : 'true'})
```

NodeJS:

```js
db = new duckdb.Database(":memory:", { allow_unsigned_extensions: "true" });
```

Secondly, you will need to set the repository endpoint in DuckDB to the HTTP url of your bucket + version of the extension
you want to install. To do this run the following SQL query in DuckDB:

```sql
SET custom_extension_repository='bucket.s3.eu-west-1.amazonaws.com/<your_extension_name>/latest';
```

Note that the `/latest` path will allow you to install the latest extension version available for your current version of
DuckDB. To specify a specific version, you can pass the version instead.

After running these steps, you can install and load your extension using the regular INSTALL/LOAD commands in DuckDB:

```sql
INSTALL duckdb_rdkit
LOAD duckdb_rdkit
```
