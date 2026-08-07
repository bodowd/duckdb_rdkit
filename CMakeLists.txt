cmake_minimum_required(VERSION 3.5)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# C++17 needs to be specified to compile rdkit
# see: https://github.com/greglandrum/rdkit_blog/blob/master/src/simple_cxx_example/CMakeLists.txt
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# Set extension name here
set(TARGET_NAME duckdb_rdkit)

set(RDKIT_PREFIX "${CMAKE_CURRENT_LIST_DIR}/build/rdkit-prefix")

find_package(RDKit QUIET)
if(NOT RDKit_FOUND)
    message(STATUS "RDKit not found; building static Boost/RDKit into ${RDKIT_PREFIX}")
    execute_process(
        COMMAND bash "${CMAKE_CURRENT_LIST_DIR}/scripts/build_rdkit_static.sh" "${RDKIT_PREFIX}"
        RESULT_VARIABLE RDKIT_BUILD_RESULT
    )
    if(NOT RDKIT_BUILD_RESULT EQUAL 0)
        message(FATAL_ERROR "scripts/build_rdkit_static.sh failed (${RDKIT_BUILD_RESULT})")
    endif()

    set(Boost_DIR "${RDKIT_PREFIX}/boost/lib/cmake/Boost-1.84.0" CACHE PATH "" FORCE)
    set(Boost_INCLUDE_DIR "${RDKIT_PREFIX}/boost/include" CACHE PATH "" FORCE)
    set(RDKit_DIR "${RDKIT_PREFIX}/rdkit/lib/cmake/rdkit" CACHE PATH "" FORCE)
    set(CMAKE_PREFIX_PATH "${RDKIT_PREFIX}/boost" CACHE PATH "" FORCE)
endif()

find_package(Boost REQUIRED COMPONENTS system serialization iostreams)
find_package(RDKit REQUIRED)

set(EXTENSION_NAME ${TARGET_NAME}_extension)
set(LOADABLE_EXTENSION_NAME ${TARGET_NAME}_loadable_extension)

project(${TARGET_NAME})
include_directories(src/include)

set(EXTENSION_SOURCES 
    src/sdf_scanner/sdf_functions.cpp
    src/sdf_scanner/sdf_scan.cpp
    src/cast.cpp
    src/mol_compare.cpp
    src/mol_formats.cpp
    src/types.cpp
    src/duckdb_rdkit_extension.cpp
    src/umbra_mol.cpp
    src/mol_descriptors.cpp
    src/qed.cpp
)

build_static_extension(${TARGET_NAME} ${EXTENSION_SOURCES})
build_loadable_extension(${TARGET_NAME} " " ${EXTENSION_SOURCES})

set(DUCKDB_RDKIT_LIBRARIES
   RDKit::SmilesParse_static
   RDKit::GraphMol_static
   RDKit::Descriptors_static
    )
# Link OpenSSL in both the static library as the loadable extension
target_link_libraries(${EXTENSION_NAME} ${DUCKDB_RDKIT_LIBRARIES})
target_link_libraries(${LOADABLE_EXTENSION_NAME} ${DUCKDB_RDKIT_LIBRARIES})


install(
  TARGETS ${EXTENSION_NAME}
  EXPORT "${DUCKDB_EXPORT_SET}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")
