FROM ubuntu:22.04
SHELL ["/bin/bash", "-c"]
RUN DEBIAN_FRONTEND=noninteractive \
  apt-get update && \ 
  sudo apt-get -y git g++ cmake ninja-build libssl-dev wget

RUN wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
  bash Miniforge3.sh -b -p "${HOME}/conda" && \ 
  source "${HOME}/conda/etc/profile.d/conda.sh" && \
  # For mamba support also run the following command
  # source "${HOME}/conda/etc/profile.d/mamba.sh" && \
  conda create -n rdkit_dev -c conda-forge -y boost-cpp boost cmake rdkit eigen librdkit-dev



#### using miniconda3 image
# FROM continuumio/miniconda3
#
#
# # install dependencies for building RDKit
# RUN conda install -c conda-forge mamba && \
#   conda create -n rdkit_dev -c conda-forge -y boost-cpp boost cmake rdkit=2023.09.4 eigen

##### building static libraries
# RUN git clone https://github.com/rdkit/rdkit.git
# RUN cd rdkit && mkdir -p build && cd build
# RUN cmake .. \
#   -DCMAKE_BUILD_TYPE=Release \
#   -DRDK_INSTALL_INTREE=ON \
#   -DRDK_BUILD_PYTHON_WRAPPERS=OFF\ 
#   -DRDK_BUILD_FUZZ_TARGETS=OFF \
#   -DRDK_INSTALL_STATIC_LIBS=ON \
#   -DBoost_USE_STATIC_LIBS=ON \
#   -DRDK_BUILD_CPP_TESTS=OFF \
#   -DBoost_NO_SYSTEM_PATHS=ON \
#   -DCMAKE_INCLUDE_PATH="${CONDA_PREFIX}/include" \
#   -DCMAKE_LIBRARY_PATH="${CONDA_PREFIX}/lib" \
#   -DCMAKE_PREFIX_PATH=$CONDA_PREFIX && \
#   make install
# try the recipe from rdkit .azure-pipelines/linux_build_fuzzer.yml




