#!/usr/bin/env python

import re
import datetime
import sys
import struct


if len(sys.argv)>1:
    file = open(sys.argv[1],"rb")


header=file.read(4)

magic=struct.unpack("4s",header)[0];

print magic;
data=file.read(12);
(n_pixels,n_emissions,n_detectors)=struct.unpack("<III",data);
print n_pixels,n_emissions,n_detectors;
while 1:
    data=file.read(8);
    if len(data)<8: 
        break;
    (lor_a,lor_b,count)=struct.unpack("HHI",data);
    print "LOR : ",(lor_a,lor_b),count;
    for item in range(count):
        data=file.read(12);
        (position,x,y,hits)=struct.unpack("IHHI",data);
        print (position,x,y,hits);
    
