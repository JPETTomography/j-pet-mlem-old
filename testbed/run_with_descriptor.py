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
args, rest = parser.parse_known_args()
recalculate = args.recalculate
print(recalculate)
print(rest)

def run_and_check(cmd):
    print("running " + ' '.join(cmd))
    info=run(cmd)
    if info.returncode !=0:
        sys.exit()

# Prepare system matrix
n_emissions = 400000
if not os.path.isfile("m_big"):
    print("m_big file does not exists: recalculating")
    recalculate=True
else:
    matrix_file = open("m_big","rb")
    matrix = SparseMatrixHeader(matrix_file)
    print(matrix.n_emissions)
    if matrix.n_emissions != n_emissions:
        recalculate=True

run(["../2d_barrel_describe_scanner","--big-barrel","-o","big_barrel"])


if recalculate :
    run_and_check(["../2d_barrel_matrix", "-c", "m_big_ref.cfg",
    "--detector-file","big_barrel_dets.txt",
    "--detector-file-sym","big_barrel_syms.txt",
    "-e", "%d" % (n_emissions,), "-o", "m_big",
         "-v"]+rest)


# Convert to full matrix
if recalculate or not os.path.isfile("f_big"):
    run_and_check(["../2d_barrel_matrix", "-c", "m_big.cfg",
    "--detector-file","big_barrel_dets.txt",
    "--detector-file-sym","big_barrel_syms.txt",
    "-o", "f_big", "-f", "m_big"])


# Prepare phantom
n_phantom_emissions = 100000000
if recalculate:
    run_and_check(["../3d_hybrid_phantom", "-c", "m_big.cfg", "-o", "p_sphere.txt",
         "-e", "%d" % (n_phantom_emissions,), "s_sphere.json", "-v"])

# Alternatively prepare phantom wih GATE 

# Reconstruct
if recalculate:
    run_and_check(["../3d_hybrid_reconstruction", "-c", "m_big.cfg", "--system", "f_big", "-o", "r_big",
         "-i", "10", "-v", "p_sphere.txt"])
