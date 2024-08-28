FROM ubuntu:22.04
SHELL ["/bin/bash", "-c"]

# install packages needed to build duckdb
RUN DEBIAN_FRONTEND=noninteractive \
  apt-get update && \ 
  sudo apt-get -y git g++ cmake ninja-build libssl-dev wget

# see: https://github.com/conda-forge/miniforge
RUN wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
  bash Miniforge3.sh -b -p "${HOME}/conda" && \ 
  source "${HOME}/conda/etc/profile.d/conda.sh" && \
  conda create -n rdkit_dev -c conda-forge -y boost-cpp boost cmake rdkit eigen librdkit-dev



