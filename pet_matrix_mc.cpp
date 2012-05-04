#include<cmath>
#include<cstdio>
#include<cstdlib>

#include<vector>

#include"event.h"
#include"detector.h"
#include"reconstruction.h"


 

const int L=32;
const int n_detectors=32;





main() {
  const int max_lors=n_detectors*n_detectors;
  double step=2.0/L;
  double *p = new double[max_lors];

  std::vector<pixel<int,double> > octant;
  for(int iy=0;iy<L/2;iy++) 
    for(int ix=iy;ix<L/2;ix++) {
      
      pixel<int, double> pix(ix,iy,step,step);
      octant.push_back(pix);
    }


  
  for(int ipix=0;ipix<octant.size();++ipix) {
    pixel<int,double> pix=octant[ipix];

    for(int i=0;i<max_lors;++i)
      p[i]=0.0;

    fill_pixel(pix,n_detectors,p,1000000);

    std::vector<Lor<int,double> > *lrs=new std::vector<Lor<int,double> >;
    double sum=0.0;
    for(int i=0;i<max_lors;++i) {
      if(p[i]>0.0) {
	int f=i/n_detectors;
	int s=i%n_detectors;
	lrs->push_back(Lor<int,double>(f,s,p[i]));
	sum+=p[i];
      }
    }

    row<pixel<int,double>, double> pixel_row(pix,*lrs);

    std::cout<<pixel_row;


  }

}


