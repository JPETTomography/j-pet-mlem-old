#ifndef __RECONSTRUCTION_H__
#define __RECONSTRUCTION_H__

#include<iostream>
#include<vector> 

#include"detector.h"

template<typename S, typename F> class Pixel {
public:
  
  typedef S short_t;

  Pixel(S i_x, S i_y, F size_x, F size_y):
    i_x_(i_x),i_y_(i_y),size_x_(size_x), size_y_(size_y) {};
					  

  S i_x() const {return i_x_;}
  S i_y() const {return i_y_;}

  F x() const {return i_x_*size_x_;}
  F y() const {return i_y_*size_y_;}

  F size_x() const {return size_x_;}
  F size_y() const {return size_y_;}

private:
  S i_x_;
  S i_y_;
 
  F size_x_;
  F size_y_;

}; 


template<typename P, typename F> class Row;

template<typename P, typename F> 
std::ostream &operator<<(std::ostream &out, const Row<P,F> &r);

template<typename P, typename F> class Row {
public:
  typedef Lor<typename P::short_t,F> lor_t;
  Row(const P &pixel, const std::vector<Lor<typename P::short_t,F> > &vec):
   pixel_(pixel), row_(vec) {}; 

   size_t fwrite(FILE *fout) const {
     size_t size_w=0;
     size_t size=row_.size();
     size_w+=::fwrite(&pixel_,sizeof(P),1,fout);
     size_w+=::fwrite(&size,sizeof(size_t),1,fout);
     size_w+=::fwrite(&row_[0],sizeof(lor_t),row_.size(),fout);
     return size_w;
   }

 private:

   size_t fread(FILE *fin) {
     size_t size_w=0;
     size_t size
     size_w+=::fread(&pixel_,sizeof(P),1,fin);
     size_w+=::fread(&size,sizeof(size_t),1,fin);
     row_.resize(size);
     size_w+=::fread(&row_[0],sizeof(lor_t),row_.size(),fout);
     return size_w;
   }

  P pixel_;
  std::vector<lor_t > row_;

   friend  std::ostream &operator<< <>(std::ostream & out, const Row &r);
  

   

};

template<typename P, typename F> 
std::ostream &operator<<(std::ostream &out, const Row<P,F> &r) {

  out<<"pixel "<<r.pixel_.i_x()<<" "<<r.pixel_.i_y();
  out<<" "<<r.pixel_.x()<<" "<<r.pixel_.y()<<" "<<r.row_.size()<<"\n";
  typename std::vector<Lor<typename P::short_t,F> >::const_iterator it=r.row_.begin();
  for(;it!=r.row_.end();++it) 
    out<<(*it).first<<" "<<(*it).second<<" "<<(*it).count<<"\n";
  out<<std::flush;
  return out;
}




template<typename F> event<F> random_from_pixel(F x, F y, F size_x, F size_y) {
  return event<F>(x+drand48()*size_x,
	       y+drand48()*size_y,
	       M_PI*drand48()
	       );
}

void fill_pixel(double x, double y, double size_x, double size_y,
		int n_detectors,
		double *p,int n) {
    for(int i=0;i<n;i++) {
      event<double> e=random_from_pixel(x,y,size_x,size_y);

    std::pair<double,double> time=tof(e,2.0);
    std::pair<short,short>   lors=lor(time,e,2.0,n_detectors);
    p[n_detectors*lors.first+lors.second]+=1.0;
    }
}

template<typename P> 
void fill_pixel(const P &pix, int n_detectors,
		double *p,int n) {
  fill_pixel(pix.x(),pix.y(),pix.size_x(),pix.size_y(),n_detectors,p,n);
}


template<typename S, typename F> class ProbabilityMatrix {
public:  
  typedef Pixel<S, F> pixel_t;
  typedef Row<pixel_t,F> row_t;

  ProbabilityMatrix(const Ring2DDetector<S,F> &scanner):scanner_(scanner) {
      for(int iy=0;iy<scanner_.n_pixels_y()/2;iy++) 
	for(int ix=iy;ix<scanner_.n_pixels_x()/2;ix++) {
      
	  octant_.push_back(
			    pixel_t(ix,iy,
				    scanner_.pixel_size_x(),
				    scanner_.pixel_size_y()
				    )
			    );
	}
  };
  
  
  S octant_size() const {return octant_.size();}
  pixel_t  octant(int i) {return octant_[i];}
  void push_back_row(row_t *row) {matrix_.push_back(row);}
  row_t *row(int i) {return matrix_[i];}

  ~ProbabilityMatrix() {
    while(!matrix_.empty()) {
      delete matrix_.back();
      matrix_.pop_back();
    }
  }


  size_t fwrite(FILE *fout) const {
    size_t size_w=0;
    for(typename std::vector<Row<pixel_t,F> *>::const_iterator it=matrix_.begin();
	it!=matrix_.end();++it) {
      size_w+=(*it)->fwrite(fout);
    }
    return size_w;
  }

 private:
  Ring2DDetector<S,F> scanner_;
  std::vector<pixel_t > octant_;
  std::vector<Row<pixel_t,F> *> matrix_;
};


template<typename S, typename F> 
  Row<Pixel<S,F>,F> *row_from_array(const Pixel<S,F> &pix, F *p,int n_detectors,
			   int max_lors) {

   std::vector<Lor<S,F> > lrs;
   lrs.reserve(4*n_detectors);
    double sum=0.0;
    for(int i=0;i<max_lors;++i) {
      if(p[i]>0.0) {
	int f=i/n_detectors;
	int s=i%n_detectors;
	lrs.push_back(Lor<S,F>(f,s,p[i]));
	sum+=p[i];
      }
    }

    return new Row<Pixel<S,F>, F>(pix,lrs);
};


#endif
