# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

### Added

- `read_sdf` to query sdf files using sql and duckdb
- sdf replacement scan to automatically detect `.sdf` and execute sql against it

## [0.2.0] - 2024-10-14

### Added

- `mol_qed` to calculate the quantitative estimate of drug-likeness (QED)

### Changed

- `mol_from_smiles` returns null if molecule cannot be made from SMILES
- Casting varchar to mol returns null if molecule cannot be made from SMILES

## [0.1.0] - 2024-08-29

### Added

- Functions to convert between SMILES and Mol objects
  - `mol_from_smiles`
  - `mol_to_smiles`
- Comparison functions
  - `is_exact_match`
  - `is_substruct`
