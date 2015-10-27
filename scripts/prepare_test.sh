#!/bin/bash
# Prepares some test input files
set -e
cd $(dirname $0)/..
set -x

# Check and change to test_input
[ -d test_input ] || mkdir test_input
cd test_input

# Scanner parameters
params='-w 0.005 -h 0.019 -d 32 -r 0.43 -p 0.01 -n 64'

# Run commands
../2d_barrel_geometry -v $params -o g_test
../2d_barrel_matrix -v $params -o m_test -e 1000000 --s-dl 0.04
../2d_barrel_matrix -c m_test.cfg -f -o f_test m_test
