#! /usr/bin/env python

from subprocess import run
import argparse

n_emissions = 10000

# Prepare system matrix
# ../2d_barrel_matrix  -c m_big.cfg  -e 1000 -v
run(["../2d_barrel_matrix" , "-c", "m_big.cfg" ,"-e", "%d" % (n_emissions,), "-o","m_big","-v"])

# Convert to full matrix
run(["../2d_barrel_matrix","-c", "m_big.cfg", "-o", "f_big", "-f","m_big"])

# Prepare phantom

# Alternatively prepare phantom wih GATE 

# Reconstruct

