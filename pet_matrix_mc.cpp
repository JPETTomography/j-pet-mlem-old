#include<cmath>
#include<cstdio>
#include<cstdlib>

#include<vector>

#include"event.h"
#include"detector.h"


template<typename F> event<F> random_from_pixel(F x, F y, F size_x, F size_y) {
  return event<F>(x+drand48()*size_x,
	       y+drand48()*size_y,
	       M_PI*drand48()
	       );
} 

const int L=128;
const int n_detectors=128;

void fill_pixel(double x, double y, double step, double *p,int n) {
    for(int i=0;i<n;i++) {
      event<double> e=random_from_pixel(x,y,step,step);

    std::pair<double,double> time=tof(e,2.0);
    std::pair<short,short> lors=lor(time,e,2.0,n_detectors);
    //    fprintf(stdout,"%4d %4d\n",lors.first, lors.second);
    p[n_detectors*lors.first+lors.second]+=1.0;
    }
}



main() {
  const int max_lors=n_detectors*n_detectors;
  double step=2.0/L;
  double *p = new double[max_lors];


  for(int iy=0;iy<L/2;iy++) 
    for(int ix=iy;ix<L/2;ix++) {
      double x=ix*step;
      double y=iy*step;
    

      for(int i=0;i<max_lors;++i)
	p[i]=0.0;

      fill_pixel(x,y,step,p,1000000);

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

      std::vector<Lor<int,double> >::iterator it=lrs->begin();
  
      fprintf(stdout,"pixel %12.8g %12.8g %5d\n",x,y,lrs->size());
      for(;it!=lrs->end();++it) {
	it->count/=sum;
	fprintf(stdout,"%4d %4d %12.8f\n",it->first,it->second,it->count);	  
      }
    }

}
