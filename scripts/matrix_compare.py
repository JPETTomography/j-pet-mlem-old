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

import petmatrix as pet

parser = argparse.ArgumentParser(description="Statistically compare two system matrices.")
parser.add_argument('--a-file', '-a', help = 'first  matrix file')
parser.add_argument('--b-file', '-b', help = 'second matrix file')
parser.add_argument('--a-emissions', type = int, help = 'number of emissions per pixel in first  matrix')
parser.add_argument('--b-emissions', type = int, help = 'number of emissions per pixel in second matrix')
parser.add_argument('--silent', '-s', action = 'store_true', help = 'prints only the significance on stdout')  
parser.add_argument('--verbose','-v', action = 'store_true') #not implemented
parser.add_argument('--graph', '-g',  action = 'store_true', help = 'displays a graph of significance for every pixels') 
args=parser.parse_args()
if args.silent:
    args.graph   = False
    args.verbose = False

if args.a_file:
    a_file = open(args.a_file,"rb")
        
    
if args.b_file:
    b_file = open(args.b_file,"rb")

a_matrix = pet.SparseMatrix(a_file)
if(args.a_emissions):
    a_matrix.header.n_emissions = args.a_emissions
a_matrix.body.Read()


b_matrix = pet.SparseMatrix(b_file)
if(args.b_emissions):
    b_matrix.header.n_emissions = args.b_emissions
b_matrix.body.Read()

a_n = a_matrix.n_emissions();
b_n = b_matrix.n_emissions();

a_p = pet.FillOctantPixMap(a_matrix.body)/a_matrix.n_emissions();
b_p = pet.FillOctantPixMap(b_matrix.body)/b_matrix.n_emissions();

diff = (a_p-b_p)
mask = (a_p+b_p)<=0

s_ab =np.sqrt((a_n*a_p*(1.0-a_p)+b_n*b_p*(1.0-b_p))/(
    (a_n + b_n -2)))*np.sqrt(1.0/a_n+1.0/b_n);
s_ab_masked =  np.ma.masked_array(s_ab,mask=mask )
n_df =  s_ab_masked.count()
tval=diff/s_ab_masked.filled(1.0);

student = stats.t(a_n+b_n-2)
p = np.ma.masked_array(2*(1.0-student.cdf(abs(tval))) , mask = mask);


chi2 = -2*np.log(p.filled(1.0)).sum();
chi2dist = stats.chi2(2*n_df)
if args.silent:
    print 1-chi2dist.cdf(chi2)
else:
    print "chi^2  = ", chi2, ", n_dof = ", 2*n_df, ", p = ", 1-chi2dist.cdf(chi2)


if args.graph:
    
    def callback(event):
        ix = math.floor(event.xdata);
        iy = math.floor(event.ydata);
        print ix, iy, pixmap[ix,iy];
    
    imgplot=plt.imshow(p.filled(0));            
    imgplot.set_interpolation("nearest")
    imgplot.figure.canvas.mpl_connect('button_press_event',callback);
    plt.colorbar(imgplot);
    plt.show()
