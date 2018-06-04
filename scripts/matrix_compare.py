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
parser.add_argument('--a-emissions', 
                    type = int, 
                    help = 'number of emissions per pixel in first  matrix'
                    )
parser.add_argument('--b-emissions', 
                    type = int, 
                    help = 'number of emissions per pixel in second matrix'
                    )
parser.add_argument('--silent', '-s', 
                    action = 'store_true', 
                    help = 'prints only the significance on stdout'
                    )  
parser.add_argument('--verbose','-v', action = 'store_true') #not implemented
parser.add_argument('--by-tors', action = 'store_true') 
parser.add_argument('--graph', '-g',  
                    action = 'store_true', 
                    help = 'displays a graph of significance for every pixels'
                    ) 

args=parser.parse_args()

if args.silent:
    args.graph   = False
    args.verbose = False

if args.a_file:
    a_file = open(args.a_file,"rb")
else: 
    print "No a file given"
    sys.exit()

if args.b_file:
    b_file = open(args.b_file,"rb")
else:
    print "No b file given"
    sys.exit()


a_matrix = pet.SparseMatrix(a_file)
if(args.a_emissions):
    a_matrix.header.n_emissions = args.a_emissions
a_matrix.body.Read()
a_matrix.body.sort_by_pixel()

b_matrix = pet.SparseMatrix(b_file)
if(args.b_emissions):
    b_matrix.header.n_emissions = args.b_emissions
b_matrix.body.Read()
b_matrix.body.sort_by_pixel()

a_n = a_matrix.n_emissions();
b_n = b_matrix.n_emissions();

print a_n
print b_n


a_p = pet.FillOctantPixMap(a_matrix.body)/a_matrix.n_emissions();
b_p = pet.FillOctantPixMap(b_matrix.body)/b_matrix.n_emissions();



diff = (a_p-b_p)
mask = (a_p+b_p)<=0

s_ab =np.sqrt(((a_n-1)*a_p*(1.0-a_p)+(b_n-1)*b_p*(1.0-b_p))/(
    (a_n + b_n -2)))*np.sqrt(1.0/a_n+1.0/b_n);
s_ab_masked =  np.ma.masked_array(s_ab,mask=mask )
n_df =  s_ab_masked.count()
tval=diff/s_ab_masked.filled(1.0);

student = stats.t(a_n+b_n-2)
p = np.ma.masked_array(2*(1.0-student.cdf(abs(tval))) , mask = mask);

#print p.filled(0);
np.savetxt("tval.txt",tval)
np.savetxt("ps.txt",p.filled(0));

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


n_d = a_matrix.header.n_detectors
n_t = a_matrix.body.n_tof_positions
n_bins = (n_d*(n_d-1))/2 * n_t
print n_bins





a_iterator = iter(a_matrix.body.items())
a_first = next(a_iterator)
b_iterator = iter(b_matrix.body.items())
b_first = next(b_iterator)
pixel_stats = []
if args.by_tors:
    while True:
        a_pixel=a_first[1]
        b_pixel=b_first[1]
        a_bins = np.zeros(n_bins)
        a_first = pet.fill_pixel_bins(n_d, n_t, a_first, a_iterator, a_bins)
        a_total = a_bins.sum()

        b_bins = np.zeros(n_bins)
        b_first = pet.fill_pixel_bins(n_d, n_t, b_first, b_iterator, b_bins)
        b_total = b_bins.sum()

        ab_bins=a_bins+b_bins;
        ab_masked_bins = np.ma.masked_equal(ab_bins, value=0)

        a_over_b = math.sqrt(float(a_total)/b_total)
        #print a_n, a_total, b_n,  b_total, a_over_b

        diff =(a_bins/a_over_b - b_bins*a_over_b)

        chi_bins = (diff*diff)/ab_masked_bins.filled(1.0);
        n_dof = ab_masked_bins.count()
        chi2 = chi_bins.sum()
        chi2dist = stats.chi2(n_dof)
        p = 1-chi2dist.cdf(chi2)
        #print a_pixel, b_pixel, " n dof = ", n_dof," chi = ", chi2, " p = ", p
        pixel_stats.append((a_pixel, n_dof, chi2,p))
        if not a_first or  not b_first:
            break
    
    ps=np.array( map(lambda item: item[3],pixel_stats)   )   
    chi2 = -2*np.log(ps).sum();
    n_dof = 2*len(ps)  
    print n_dof, chi2, 1-stats.chi2.cdf(chi2, n_dof)        
                            
                            
            
