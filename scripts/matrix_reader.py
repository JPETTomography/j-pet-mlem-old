#!/usr/bin/env python

import re
import datetime
import sys

import math

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np



import petmatrix as pet


if len(sys.argv)>1:
    file = open(sys.argv[1],"rb")

matrix = pet.SparseMatrix(file)
matrix.show();
matrix.body.Read()
print matrix.body.stats()[0:3]


def callback(event):
    ix = math.floor(event.xdata);
    iy = math.floor(event.ydata);
    print ix, iy, pixmap[ix,iy];
    
    


pixmap = pet.FillPixMap(matrix.body)
imgplot=plt.imshow(pixmap)
imgplot.set_interpolation("nearest")
imgplot.figure.canvas.mpl_connect('button_press_event',callback);
plt.colorbar(imgplot);
plt.show()
