#!/usr/bin/env bash
set -euo pipefail

BOOST_VERSION="1.84.0"
BOOST_GIT_TAG="boost-1.84.0"
BOOST_GIT_URL="https://github.com/boostorg/boost.git"

RDKIT_TAG="Release_2024_03_5"
RDKIT_VERSION="2024.03.5"
RDKIT_TARBALL="${RDKIT_TAG}.tar.gz"
RDKIT_URL="https://github.com/rdkit/rdkit/archive/refs/tags/${RDKIT_TARBALL}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${1:-${REPO_ROOT}/build/rdkit-prefix}"

SRC_DIR="${PREFIX}/src"
BOOST_SRC_DIR="${SRC_DIR}/boost-${BOOST_VERSION}"
BOOST_PREFIX="${PREFIX}/boost"
RDKIT_PREFIX="${PREFIX}/rdkit"
BOOST_BUILD_DIR="${SRC_DIR}/boost-build"
RDKIT_BUILD_DIR="${SRC_DIR}/rdkit-build"

find_rdkit_cmake_dir() {
    find "${RDKIT_PREFIX}" -type f -name rdkit-config.cmake -path '*/cmake/rdkit/*' 2>/dev/null | head -1 | xargs dirname 2>/dev/null
}

if [ -n "$(find_rdkit_cmake_dir)" ]; then
    echo "Static RDKit already installed at ${RDKIT_PREFIX}; skipping."
    exit 0
fi

for tool in cmake curl tar find git; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "error: required tool '${tool}' not found" >&2
        exit 1
    fi
done

if command -v nproc >/dev/null 2>&1; then
    JOBS="$(nproc)"
elif command -v sysctl >/dev/null 2>&1; then
    JOBS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
else
    JOBS=4
fi

download_and_extract() {
    local url="$1" tarball="$2" outdir="$3"
    if [ ! -f "${SRC_DIR}/${tarball}" ]; then
        echo "Downloading ${url}"
        curl -L --fail --retry 3 -o "${SRC_DIR}/${tarball}" "${url}"
    fi
    if [ ! -d "${outdir}" ]; then
        echo "Extracting ${tarball}"
        tar -xf "${SRC_DIR}/${tarball}" -C "${SRC_DIR}"
    fi
}

mkdir -p "${SRC_DIR}"

echo "==> Cloning Boost ${BOOST_VERSION} (static, PIC)"
if [ ! -d "${BOOST_SRC_DIR}" ]; then
    git clone --depth 1 --branch "${BOOST_GIT_TAG}" \
        --recurse-submodules --shallow-submodules \
        "${BOOST_GIT_URL}" "${BOOST_SRC_DIR}"
fi
cmake -S "${BOOST_SRC_DIR}" -B "${BOOST_BUILD_DIR}" \
    -U BOOST_INCLUDE_LIBRARIES \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${BOOST_PREFIX}" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
cmake --build "${BOOST_BUILD_DIR}" --target install -j "${JOBS}"

echo "==> Building RDKit ${RDKIT_VERSION} (static, PIC)"
download_and_extract "${RDKIT_URL}" "${RDKIT_TARBALL}" "${SRC_DIR}/rdkit-${RDKIT_TAG}"
cmake -S "${SRC_DIR}/rdkit-${RDKIT_TAG}" -B "${RDKIT_BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${RDKIT_PREFIX}" \
    -DCMAKE_PREFIX_PATH="${BOOST_PREFIX}" \
    -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
    -DRDK_BUILD_CPP_TESTS=OFF \
    -DRDK_BUILD_DESCRIPTORS3D=OFF \
    -DRDK_BUILD_FREETYPE_SUPPORT=OFF \
    -DRDK_INSTALL_STATIC_LIBS=ON \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DCMAKE_C_STANDARD=17 \
    -DRDK_INSTALL_INTREE=OFF \
    -DBUILD_SHARED_LIBS=OFF
cmake --build "${RDKIT_BUILD_DIR}" -j "${JOBS}"
cmake --install "${RDKIT_BUILD_DIR}"

RDKIT_CMAKE_DIR="$(find_rdkit_cmake_dir)"
if [ -z "${RDKIT_CMAKE_DIR}" ]; then
    echo "error: build finished but rdkit-config.cmake was not installed" >&2
    exit 1
fi

BOOST_CMAKE_DIR="$(find "${BOOST_PREFIX}" -type f -name BoostConfig.cmake 2>/dev/null | head -1 | xargs dirname 2>/dev/null)"

echo
echo "Done. Static Boost/RDKit installed under ${PREFIX}."
echo "Pass these EXT_FLAGS to the duckdb_rdkit build:"
echo "  make release GEN=ninja EXT_FLAGS=\"-DRDKit_DIR=${RDKIT_CMAKE_DIR} \\"
echo "           -DBoost_DIR=${BOOST_CMAKE_DIR} \\"
echo "           -DBoost_INCLUDE_DIR=${BOOST_PREFIX}/include \\"
echo "           -DCMAKE_PREFIX_PATH=${BOOST_PREFIX} -U boost_assert_DIR\""
