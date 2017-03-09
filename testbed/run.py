#! /usr/bin/env python

import os.path
import sys

sys.path.append("../scripts")

from petmatrix import SparseMatrixHeader

from subprocess import run
import argparse

recalculate = False

parser = argparse.ArgumentParser(description="Full reconstruction workflow")
parser.add_argument('--recalculate', '-r', action='store_true', dest='recalculate')
args = parser.parse_args()
recalculate = args.recalculate


# Prepare system matrix
n_emissions = 1000000
if not os.path.isfile("m_bin"):
    recalculate=True
else:
    matrix = SparseMatrixHeader("m_bin")
    print(matrix.n_emissions)

if recalculate:
    run(["../2d_barrel_matrix", "-c", "m_big.cfg", "-e", "%d" % (n_emissions,), "-o", "m_big",
         "-v"])

# Convert to full matrix
run(["../2d_barrel_matrix", "-c", "m_big.cfg", "-o", "f_big", "-f", "m_big"])

# Prepare phantom
n_phantom_emissions = 100000000
if recalculate:
    run(["../3d_hybrid_phantom", "-c", "m_big.cfg", "-o", "p_sphere.txt",
         "-e", "%d" % (n_phantom_emissions,), "s_sphere.json", "-v"])

# Alternatively prepare phantom wih GATE 

# Reconstruct
if recalculate:
    run(["../3d_hybrid_reconstruction", "-c", "m_big.cfg", "--system", "f_big", "-o", "r_big",
         "-i", "10", "-v", "p_sphere.txt"])
