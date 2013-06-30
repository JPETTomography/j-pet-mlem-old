#pragma once



template<typename F>
class Slice {
public:
  Slice(F *rep, int offset, int stride) 
    : rep_(rep), 
      offset_(offset), 
      stride_(stride) {}

  F &operator[](int i) {return rep_[offset_+stride_*i];}
  F  operator[](int i) const {return rep_[offset_+stride_*i];}

private:
  F *rep_;
  int offset_;
  int stride_;
};


template<int N, typename F>
class Vector {
public:
  Vector() {};
  Vector(F *rep) {
    for(int i=0;i<N;i++) 
      rep_[i]=rep[i];
  }

  template<typename I> Vector(I begin) {
    I it=begin;
    for(int i=0;i<N;++i,++it) 
      rep_[i]=*it;
  }


  F &operator[](int i) {return rep_[i];};
  F  operator[](int i) const {return rep_[i];};
  
private:
  F rep_[N];
};


template<int N, int M, typename F = double > 
class Matrix {
public:
  Matrix() {};
  Matrix(F *rep) {
    for(int i=0;i<N*M;++i) {
      rep_[i] = rep[i];
    }
  }
  template<typename I> Matrix(I begin) {
    I it=begin;
    for(int i=0;i<N*M;++i,++it) {
	rep_[i] = *it;
      }
  };
  
  Slice<F> operator[](int i) {
    return Slice<F>(rep_,i*M,1);
  }
private:
  F rep_[N*M];
};
