#ifndef __PIXEL_GRID__
#define __PIXEL_GRID__

#include<vector>

template<typename F=double> class Point {
 public:
  Point(F xa, F ya):x(xa),y(ya){}    
  F x;
  F y;
};

template<typename I=int> class Index {
 public:
  Index(I xa, I ya):x(xa),y(ya) {}    
  I x;
  I y;
};

template<typename F = double, typename I = int >
class PixelGrid {
public:

  PixelGrid(const Point<F> &ll, const Point<F> &ur, I nx, I ny):
    nx_(nx),ny_(ny),
    pixels_(nx_*ny_,(F)0.0) {
    ll_x_=ll.x;
    ll_y_=ll.y;
    ur_x_=ur.x;
    ur_y_=ur.y;

    dx_=(ur_x_-ll_x_)/nx_;
    dy_=(ur_y_-ll_y_)/ny_;
    
  };


  I nx() const {return nx_;}
  I ny() const {return ny_;}

  Point<F> ll() const {return Point<F>(ll_x_,ll_y_);}
  Point<F> ur() const {return Point<F>(ur_x_,ur_y_);}
  

  F dx() const {return dx_;}
  F dy() const {return dy_;}


  Point<F> center(I ix, I iy) const  {
    return Point<F>(ll_x_+(ix+0.5)*dx_,ll_y_+(iy+0.5)*dy_);
  }
  Point<F> center(const Index<I> &ind) const  {
    return center(ind.x,ind.y);
  }

  Index<I> in(F x, F y) const {
    I ix=(int) floor( (x-ll_x_)/dx_);
    I iy=(int) floor( (y-ll_y_)/dy_);
    return Index<I>(ix,iy);
  };

  Index<I> in(const Point<F> &p) const {
    return in(p.x,p.y);
  };
  
  void add(I ix, I iy, F w=1.0) {
    pixels_[index(ix,iy)]+=w;
  }
  
  void add(Index<I> ind, F w=1.0) {
    add(ind.x,ind.y,w);
  }

  void insert(F x, F y, F w) {
    add(in(x,y),w);
  }

  void insert(const Point<F> &p, F w) {
    add(in(p),w);
  }

  I index(I ix, I iy) const {return iy*nx_+ix;}
  I index(Index<F> ind) const {return index(ind.x,ind.y);}


  F operator()(I i) const {return pixels_[i];}
  F operator[](I i) const {return this->operator()(i);}
  F operator()(I ix,I iy) const {return pixels_[index(ix,iy)];}

private:
  F  ll_x_, ll_y_;
  F  ur_x_, ur_y_;
  I  nx_, ny_;
  
  std::vector<F> pixels_;
  F dx_,dy_;
  
  

};

#endif
