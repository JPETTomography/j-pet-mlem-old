import struct
import numpy as np

def index(row, column):
    return row*n_pixels+column;
def key(elem):
    return index(elem[1][1],elem[1][2]);


class SparseMatrixHeader(object):
    def __init__ (self, file):
        self.magick = struct.unpack("4s",file.read(4))[0]
        data=file.read(12);
        (self.n_pixels,self.n_emissions,self.n_detectors) = struct.unpack("<III",data)
        
    def show(self):
        print (self.magick)
        print (self.n_pixels, self.n_detectors, self.n_emissions)


class SparseMatrixBody(object):

    def __init__(self):
        self._items=[]
    
    def Read(self):
        while self.ReadTOR() >0 :
            pass
                
        
    def stats(self):
        counter = 0
        total = 0
        max = 0
        for item in self._items:
            counter = counter + 1
            total = total + item[2]
            if max < item[2]:
                max = item[2];
                
            hist = [0]*(max+1)
            
        for item in self._items:
            hist[item[2]]=hist[item[2]]+1
                                                   
        return (counter, total, max, hist)
       
    def items(self):
        return self._items

    def sort_by_pixel(self):
        self._items.sort(key=lambda item: item[1])
     
        
    def show(self):
        print (self.n_tof_positions)

class SparseMatrixTOFpBody(SparseMatrixBody):
    def __init__(self, header, file):
        super(SparseMatrixTOFpBody,self).__init__()
        self.header = header;
        self.header.n_pixels = 2*self.header.n_pixels
        self.file = file
        self.n_tof_positions = struct.unpack("<I",file.read(4))[0]           
       
            
            
    def ReadTOR(self):
        TOR = []
        data = self.file.read(4)
        if len(data) < 4:
            return 0
        (a,b) =  struct.unpack("<HH",data) 
        count =  struct.unpack("<I",self.file.read(4))[0]
        
        for i in range(count):
            data = self.file.read(12);
            hit = struct.unpack("IHHI",data) #(tof, x, y , count)    
            self._items.append( ( (a,b,hit[0]), (hit[1], hit[2]),hit[3])  )
        return count;
            
class SparseMatrixPETpBody(SparseMatrixBody):
    def __init__(self, header, file):
        super(SparseMatrixPETpBody,self).__init__()
        self.header = header;
        self.header.n_pixels = 2*self.header.n_pixels
        self.file = file
        self.n_tof_positions = 1           
      
            
            
    def ReadTOR(self):
        TOR = []
        data = self.file.read(4)
        if len(data) < 4:
            return 0
        (a,b) =  struct.unpack("<HH",data) 
        count =  struct.unpack("<I",self.file.read(4))[0]        
        for i in range(count):
            data = self.file.read(8);
            hit = struct.unpack("HHI",data) #(tof, x, y , count)    
            self._items.append( ( (a,b,0), (hit[0], hit[1]),hit[2])  )
        return count;        
    
            
class SparseMatrix:
    def __init__(self,file):
        self.header = SparseMatrixHeader(file);
        if self.header.magick == "TOFp":
            self.body = SparseMatrixTOFpBody(self.header, file)
        elif self.header.magick == "PETp":
            self.body = SparseMatrixPETpBody(self.header, file)
    
    def n_emissions(self):
        return self.header.n_emissions
            
    def show(self):
        self.header.show();
        self.body.show();

            
def FillOctantPixMap(matrix):
    pixmap=np.zeros((matrix.header.n_pixels/2,matrix.header.n_pixels/2));

    for item in matrix.items():               
        (ix,iy)  = item[1]
        count  = item[2] #if ix!= iy else 2*item[2]
        pixmap[ix,iy]+=count
    return pixmap    
        
def FillPixMap(matrix):
    pixmap=np.zeros((matrix.header.n_pixels,matrix.header.n_pixels));
    center = matrix.header.n_pixels/2;
    for item in matrix.items():
               
        (x,y)  = item[1]
        count  = item[2]

        ix = x + center;     
               
        iy = y + center;               
        pixmap[ix,iy]+=count
        pixmap[iy,ix]+=count                    
               
        iy = center -  y -1;                    
        pixmap[ix,iy]+=count
        pixmap[iy,ix]+=count
                
        ix = center - x - 1;
                     
        iy = y+center;                    
        pixmap[ix,iy]+=count
        pixmap[iy,ix]+=count  
                  
        iy = center -  y - 1;                    
        pixmap[ix,iy]+=count
        pixmap[iy,ix]+=count
                    
    return pixmap        




def offset_in_tor(n_detectors, n_tof_pos, tor):
    (d1, d2, dl) = tor
    offset = (d1*(d1-1))/2
    return (offset + d2)*n_tof_pos +  dl

def bin_item(n_detectors, n_tof_pos, item, bins):
    (tor, pixel, hits) =  item
    index = offset_in_tor(n_detectors, n_tof_pos, tor)
    bins[index] += hits
    #print index, hits

def fill_pixel_bins(n_detectors, n_tof_pos, first_item, iter_, bins):
    bin_item(n_detectors, n_tof_pos, first_item, bins) 
    pixel = first_item[1]
    while True:
        try:
            item = next(iter_)
        except StopIteration:
            return False
        if(item[1] != first_item[1]):
            return item
        bin_item(n_detectors, n_tof_pos, item, bins)
        
        

    
