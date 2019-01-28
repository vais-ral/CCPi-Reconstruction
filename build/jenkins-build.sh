#!/usr/bin/env bash
export CCPI_PREBUILD="recipes/library -c conda-forge -c ccpi"
bash <(curl -L https://raw.githubusercontent.com/vais-ral/CCPi-VirtualMachine/master/scripts/jenkins-build.sh)
