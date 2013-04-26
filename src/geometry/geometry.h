#pragma once

namespace Geometry {

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
