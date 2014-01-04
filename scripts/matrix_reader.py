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


class SparseMatrixHeader:
    def __init__ (self, file):
        self.magick = struct.unpack("4s",file.read(4))[0]
        data=file.read(12);
        (self.n_pixels,self.n_emissions,self.n_detectors) = struct.unpack("<III",data)
        
    def show(self):
        print self.magick
        print self.n_pixels, self.n_detectors, self.n_emissions    


class SparseMatrixTOFpBody:
        def __init__(self, header, file):
            self.header = header;
            self.header.n_pixels = 2*self.header.n_pixels
            self.file = file
            self.n_tof_positions = struct.unpack("<I",file.read(4))[0]
            self.TORS=[];
            
            
        def ReadTOR(self):
            TOR = []
            data = file.read(4)
            if len(data) < 4:
                return []
            (a,b) =  struct.unpack("<HH",data) 
            count =  struct.unpack("<I",file.read(4))[0]
            #print a,b,count
            for i in range(count):
                data = file.read(12);
                hit = struct.unpack("IHHI",data);
                #print hit
                TOR.append(hit)
            return TOR;
            
        def Read(self):
            while 1:
                tor = self.ReadTOR()
                if tor == []:
                    break
                self.TORS.append(tor)    
        
        def FillPixMap(self):
            pixmap=np.zeros((self.header.n_pixels,self.header.n_pixels));
            center = self.header.n_pixels/2;
            for tor in self.TORS:
                for hit in tor:
                    ix = hit[1]+center;
                    
                    iy = hit[2]+center;                    
                    pixmap[ix,iy]+=hit[3]
                    pixmap[iy,ix]+=hit[3]                    
                    iy = center -  hit[2]-1;                    
                    pixmap[ix,iy]+=hit[3]
                    pixmap[iy,ix]+=hit[3]
                    
                    ix = center - hit[1] - 1;
                     
                    iy = hit[2]+center;                    
                    pixmap[ix,iy]+=hit[3]
                    pixmap[iy,ix]+=hit[3]                    
                    iy = center -  hit[2]-1;                    
                    pixmap[ix,iy]+=hit[3]
                    pixmap[iy,ix]+=hit[3]
                    
            return pixmap        
        
        def stats(self):
            counter = 0
            total = 0
            max = 0
            for tor in self.TORS:
                for hit in tor:
                    counter = counter + 1
                    total = total + hit[3]
                    if max < hit[3]:
                        max = hit[3];  
            hist = [0]*(max+1)
            for tor in self.TORS:
                for hit in tor:
                    hist[hit[3]]=hist[hit[3]]+1
                                                   
            return (counter, total, max, hist)
            
        def show(self):
            print self.n_tof_positions
    
            
class SparseMatrix:
    def __init__(self,file):
        self.header = SparseMatrixHeader(file);
        if self.header.magick == "TOFp":
            self.body = SparseMatrixTOFpBody(self.header, file)
        
    def show(self):
        self.header.show();
        self.body.show();
            

if len(sys.argv)>1:
    file = open(sys.argv[1],"rb")

matrix = SparseMatrix(file)
matrix.show();
matrix.body.Read()
print matrix.body.stats()

pixmap=matrix.body.FillPixMap()
imgplot=plt.imshow(pixmap)
imgplot.set_interpolation("nearest")
plt.colorbar(imgplot);
plt.show()