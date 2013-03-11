#!/usr/bin/env python

import re
import datetime
import sys
import struct

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np


def index(row, column):
    return row*n_pixels+column;
def key(elem):
    return index(elem[1][1],elem[1][2]);

if len(sys.argv)>1:
    file = open(sys.argv[1],"rb")


header=file.read(4)

magic=struct.unpack("4s",header)[0];

print magic;
data=file.read(12);
(n_pixels,n_emissions,n_detectors)=struct.unpack("<III",data);
print n_pixels,n_emissions,n_detectors;
matrix = [];
while 1:
    data=file.read(8);
    if len(data)<8:
        break;
    (lor_a,lor_b,count)=struct.unpack("HHI",data);
#    print "LOR : ",(lor_a,lor_b),count;
    for item in range(count):
        data=file.read(12);
        (position,x,y,hits)=struct.unpack("IHHI",data);
#        print (position,x,y,hits);
        matrix.append(((lor_a,lor_b),(position,x,y,hits)));

matrix.sort(key=key);
print len(matrix);





pixmap=np.zeros((n_pixels,n_pixels));

for element in matrix:
    x=element[1][1];
    y=element[1][2];
    pixmap[x,y]+=element[1][3];


imgplot=plt.imshow(pixmap)
imgplot.set_interpolation("nearest")
plt.show()
