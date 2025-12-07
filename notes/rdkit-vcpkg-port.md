# RDKit vcpkg Port Notes

## Summary

RDKit is not available in vcpkg. Creating a C++ only port is feasible.

## Current Build Issue

The project uses:
- **vcpkg toolchain** for DuckDB extension build system
- **Spack-installed RDKit** at `/mnt/aux-data/teague/Dev/spack/var/spack/environments/duckdb/.spack-env/view/`

The conflict: RDKit's cmake config (`rdkit-targets.cmake:61-64`) declares dependencies on `Boost::system`, `Boost::serialization`, `Boost::iostreams`. When vcpkg's toolchain intercepts `find_package(Boost)`, it looks in vcpkg's install tree where Boost doesn't exist.

## RDKit Dependencies (from Spack recipe)

### Required
- `boost` (+system +serialization +iostreams) - **in vcpkg**
- `sqlite` - **in vcpkg**

### Optional (for full build)
- `freetype` - **in vcpkg**
- `eigen3` (for 3D descriptors) - **in vcpkg**
- `coordgen` - **NOT in vcpkg** (would need port)
- `maeparser` - **NOT in vcpkg** (would need port)
- `freesasa` - **NOT in vcpkg**
- Python/NumPy (for wrappers) - not needed for C++ only

### For C++ Only Build
Only need: boost, sqlite, optionally eigen3/freetype. All available in vcpkg.

## vcpkg Port Structure

A port requires two files in `/home/teague/Dev/vcpkg/ports/rdkit/`:

### vcpkg.json
```json
{
  "name": "rdkit",
  "version": "2024.03.3",
  "description": "RDKit: Open-Source Cheminformatics Software",
  "homepage": "https://www.rdkit.org",
  "license": "BSD-3-Clause",
  "dependencies": [
    "boost-system",
    "boost-serialization",
    "boost-iostreams",
    "sqlite3",
    {
      "name": "vcpkg-cmake",
      "host": true
    },
    {
      "name": "vcpkg-cmake-config",
      "host": true
    }
  ],
  "features": {
    "freetype": {
      "description": "Build with FreeType support",
      "dependencies": ["freetype"]
    },
    "3d": {
      "description": "Build 3D descriptor calculators",
      "dependencies": ["eigen3"]
    }
  }
}
```

### portfile.cmake (skeleton)
```cmake
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO rdkit/rdkit
    REF Release_2024_03_3
    SHA512 <calculate-sha512>
    HEAD_REF master
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DRDK_INSTALL_INTREE=OFF
        -DRDK_BUILD_PYTHON_WRAPPERS=OFF
        -DRDK_BUILD_COORDGEN_SUPPORT=OFF
        -DRDK_BUILD_MAEPARSER_SUPPORT=OFF
        -DRDK_BUILD_FREESASA_SUPPORT=OFF
        -DRDK_BUILD_YAEHMOP_SUPPORT=OFF
        -DRDK_BUILD_XYZ2MOL_SUPPORT=OFF
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(CONFIG_PATH lib/cmake/rdkit)
vcpkg_copy_pdbs()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/license.txt")
```

## Effort Estimate (C++ Only)

| Task | Time |
|------|------|
| Create vcpkg.json | 15 min |
| Create portfile.cmake | 30 min |
| Debug build issues | 1-2 hours |
| Test on Linux | 30 min |
| **Total** | **2-4 hours** |

## Alternative: Fix Current Build

Instead of creating a port, could:
1. Add boost to vcpkg manifest so `Boost::system` target exists
2. Disable vcpkg toolchain's find_package override for specific packages
3. Use Spack entirely (remove vcpkg toolchain)

## References

- Spack recipe: `spack edit rdkit`
- RDKit CMake docs: https://github.com/rdkit/rdkit/blob/master/Code/cmake/Modules/
- vcpkg port tutorial: https://learn.microsoft.com/en-us/vcpkg/get_started/get-started-adding-to-registry
- Example complex port: `/home/teague/Dev/vcpkg/ports/eigen3/`
