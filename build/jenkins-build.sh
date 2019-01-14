#!/usr/bin/env bash

if [[ -n ${CIL_VERSION} ]]
then
  echo Using defined version: $CIL_VERSION
else
  export CIL_VERSION=0.10.4
  echo Defining version: $CIL_VERSION
fi
# Script to builds source code in Jenkins environment
# module try-load conda

# install miniconda if the module is not present
if hash conda 2>/dev/null; then
  echo using conda
else
  if [ ! -f Miniconda3-latest-Linux-x86_64.sh ]; then
    wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
  fi
  ./Miniconda3-latest-Linux-x86_64.sh -u -b -p .
  PATH=$PATH:./bin
fi

# presume that git clone is done before this script is launched, if not, uncomment
# git clone https://github.com/vais-ral/CCPi-Regularisation-Toolkit
conda install -y conda-build
#export CIL_VERSION=0.10.2
#cd CCPi-Regularisation-Toolkit # already there by jenkins
# Reconstruction

conda build recipes/library -c conda-forge -c ccpi
export REG1_FILES=`conda build Wrappers/Python/conda-recipe --output`
# REG1_FILES variable contains output files of this build

# need to call first build
conda build Wrappers/Python/conda-recipe
# then need to call the same with --output 
export REG_FILES=`conda build Wrappers/Python/conda-recipe --output`
# REG_FILES variable contains output files
# upload files to anaconda
if [[ -n ${CCPI_CONDA_TOKEN} ]]
then
  conda install anaconda-client
  while read -r outfile; do
  anaconda -v -t ${CCPI_CONDA_TOKEN}  upload $outfile --force --label dev
  done <<< "$REG1_FILES"
  while read -r outfile; do
  anaconda -v -t ${CCPI_CONDA_TOKEN}  upload $outfile --force --label dev
  done <<< "$REG_FILES"
else
  echo CCPI_CONDA_TOKEN not defined, will not upload to anaconda.
fi
