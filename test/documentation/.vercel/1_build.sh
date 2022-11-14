#!/usr/bin/env bash
set -euxo pipefail

WORK_DIR=`pwd`
EXPORT_DIR="${WORK_DIR}/export"      # ${EXPORT_DIR}/html will be served by vercel.
PATH="${WORK_DIR}/doxygen-bin:$PATH" # Add binaries from 0_setup.sh to PATH.
CACHE_DIR="${WORK_DIR}/node_modules" # The node_modules directory is always cached.
BUILD_DIR="${CACHE_DIR}/build"       # We also want to cache the documentation build.

### Print versions.
cmake3 --version
doxygen --version

### Configure documentation build.
mkdir -p "${BUILD_DIR}" && cd "${BUILD_DIR}"
cmake3 "${WORK_DIR}/.." \
    -DSEQAN3_VERCEL_PREVIEW_DOC=ON \
    -DSEQAN3_VERCEL_URL="${VERCEL_URL}" \
    -DCMAKE_INSTALL_PREFIX="" \
    -DCMAKE_INSTALL_DOCDIR="."  1>/dev/null

### Build documentation.
cmake3 --build . --target download-cppreference-doxygen-web-tag 1>/dev/null
ctest3 -j4 --output-on-failure .

### Install documentation.
mkdir -p "${EXPORT_DIR}/usr/" "${EXPORT_DIR}/dev/"
DESTDIR="${EXPORT_DIR}/usr/" cmake3 -DCOMPONENT=doc -P cmake_install.cmake 1>/dev/null
DESTDIR="${EXPORT_DIR}/dev/" cmake3 -DCOMPONENT=doc-dev -P cmake_install.cmake 1>/dev/null
cp "${SOURCE_DIR}/test/documentation/.vercel/index.html" "${EXPORT_DIR}/index.html"

### Run indexer.
# We want the resulting index to be in the binary dir for easy access in api/doxysearch.sh.
mkdir -p "${WORK_DIR}/doxygen-bin/usr" "${WORK_DIR}/doxygen-bin/dev"
# Will put a directory doxysearch.db in the current directory.
cd ${WORK_DIR}/doxygen-bin/usr
doxyindexer "${BUILD_DIR}/usr/searchdata.xml"
cd ${WORK_DIR}/doxygen-bin/dev
doxyindexer "${BUILD_DIR}/dev/searchdata.xml"
