#pragma once

namespace Geometry {

  namespace Base {
  template<int N, typename F = double> class Vector {
  public:
    Vector(F v[N]){ 
      for(int i=0;i<N;i++) {
        v_[i]=v[i];
      }
    }
    
    Vector &operator+=(const Vector &rhs) {
      for(int i=0;i<N;++i)
        v_[i]+=rhs[i];

      return *this;
    };
    Vector &operator-=(const Vector &rhs) {
      for(int i=0;i<N;++i)
        v_[i]-=rhs[i];
      return *this;
    };

    Vector &operator+=(F rhs) {
      for(int i=0;i<N;++i)
        v_[i]+=rhs;
      return *this;
    };
    Vector &operator-=(F rhs) {
      for(int i=0;i<N;++i)
        v_[i]-=rhs;
      return *this;
    };

    Vector &operator*=(F rhs){
      for(int i=0;i<N;++i)
        v_[i]*=rhs;
      return *this;
    };;
    Vector &operator/=(F rhs){
      for(int i=0;i<N;++i)
        v_[i]/=rhs;
      return *this;
    };

    
    F &operator[](int i) {return v_[i];}
    F  operator[](int i) const  {return v_[i];}
    

  private:
    F v_[N];
  };

  template<int N, typename F> class Point { 
  public:
    Point(F p[N]) {
      for(int i=0;i<N;++i)
          p_[i]+=p[i];        
    }


    Point &operator+=(const Vector<N,F> &rhs) {
      for(int i=0;i<N;++i)
        p_[i]+=rhs[i];

      return *this;
    };
    Point &operator-=(const Vector<N,F> &rhs) {
      for(int i=0;i<N;++i)
        p_[i]-=rhs[i];
      return *this;
    };
    
  private:
    F p_[N];
  };



  template<int N, typename F=double> 
  Point<N,F> operator+(const Vector<N,F> &v, const Point<N,F> &p);
  template<int N, typename F=double> 
  Point<N,F> operator-(const Vector<N,F> &v, const Point<N,F> &p);
  template<int N, typename F=double> 
  Point<N,F> operator+( const Point<N,F> &p, const Vector<N,F> &v);  
  template<int N, typename F=double> 
  Point<N,F> operator+( const Point<N,F> &p, const Vector<N,F> &v);

  template<int N, typename F=double>  
  Vector<N,F> operator+(const Vector<N,F> &v1, const Vector<N,F> &v2);
  template<int N, typename F=double> 
  Vector<N,F> operator-(const Vector<N,F> &v1, const Vector<N,F> &v2);

  template<int N, typename F=double> 
  Vector<N,F> operator-(const Point<N,F> &p1, const Point<N,F> &p2);
  

  template<int N, typename F=double> 
  Vector<N,F> operator+(const Vector<N,F> &v1, F s);
  template<int N, typename F=double> 
  Vector<N,F> operator+(F s, const Vector<N,F> &v1);
  }


  template<int N, typename F=double> class Vector: public Base::Vector<N,F> {
  public:
    Vector(F v[N]):Base::Vector<N,F>(v) {};
  };
  template<int N, typename F=double> class Point: public Base::Point<N,F> {
  public:
    Point(F v[N]):Base::Point<N,F>(v) {};
  };

  template<typename F> class Vector<2,F>: public Base::Vector<2,F> {
  public:
    Vector(F v[2]):Base::Vector<2,F>(v) {};

    Vector &rotate(F angle) {
      F s=sin(angle);
      F c=cos(angle);
      F x=(*this)[0];
      F y=(*this)[1];

      (*this)[0]=x *c - y *s;
      (*this)[1]=y *c + x *s;
      
    }
  };

  template<typename F> class Point<2,F>: public Base::Point<2,F> {
  public:
    Point(F v[2]):Base::Point<2,F>(v) {};

    Point &rotate(F angle, Point center) {
      Vector<2,F> v=(*this)-center;
      v.rotate(angle);
      Point p=center+v;
      (*this)[0]=p[0];
      (*this)[1]=p[1];

      return *this;
    }
    Point &rotate(F angle) {
      Vector<2,F> v;
      v[0]=(*this)[0];
      v[1]=(*this)[1];
      v.rotate(angle);
      (*this)[0]=v[0];
      (*this)[1]=v[1];

      return *this;
    }
  };



  

}
