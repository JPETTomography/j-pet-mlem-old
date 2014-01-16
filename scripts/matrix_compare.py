#!/usr/bin/env python

import re
import datetime
import sys

import math

import argparse

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import scipy as sp
import scipy.stats as stats
from scipy.stats import t


import petmatrix as pet

parser = argparse.ArgumentParser(description='Statistically compare two system matrices.')
parser.add_argument('--a-file', '-a')
parser.add_argument('--b-file', '-b')
args=parser.parse_args()


if args.a_file:
    a_file = open(args.a_file,"rb")
if args.b_file:
    b_file = open(args.b_file,"rb")

a_matrix = pet.SparseMatrix(a_file)
a_matrix.show();
a_matrix.body.Read()
print a_matrix.body.stats()[0:2]
print a_matrix.n_emissions()

b_matrix = pet.SparseMatrix(b_file)
b_matrix.show();
b_matrix.body.Read()
print b_matrix.body.stats()[0:2]
print b_matrix.n_emissions()

def callback(event):
    ix = math.floor(event.xdata);
    iy = math.floor(event.ydata);
    print ix, iy, pixmap[ix,iy];
    
    

a_n=a_matrix.n_emissions();
b_n=b_matrix.n_emissions();

a_p = pet.FillPixMap(a_matrix.body)/a_matrix.n_emissions();
b_p = pet.FillPixMap(b_matrix.body)/b_matrix.n_emissions();
diff = (a_p-b_p);
s_ab =np.sqrt((a_n*a_p*(1.0-a_p)+b_n*b_p*(1.0-b_p))/(
    (a_matrix.n_emissions() +     a_matrix.n_emissions() -2)))*np.sqrt(1.0/a_n+1.0/b_n);

imgplot=plt.imshow(diff/s_ab)

imgplot.set_interpolation("nearest")
imgplot.figure.canvas.mpl_connect('button_press_event',callback);
plt.colorbar(imgplot);
plt.show()
