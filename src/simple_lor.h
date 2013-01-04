#pragma once
#include<iostream>
#include<algorithm>

class SimpleLor {
public:

 SimpleLor(int first_a = 0,int second_a = 0):first(first_a),second(second_a){};


  int first;
  int second;

  static int t_index(const SimpleLor &lor) {
    SimpleLor ordered_lor=lor;
    if(ordered_lor.first< ordered_lor.second)
      std::swap(ordered_lor.first, ordered_lor.second);
    return (ordered_lor.first*(ordered_lor.first+1))/2+ordered_lor.second;
  }

  static void init(int n_detectors_a) {
    n_detectors_=n_detectors_a;
    n_lors_     =(n_detectors_*(n_detectors_+1))/2;
  }

  static int n_lors() {return n_lors_;}

  bool operator==(const SimpleLor &rhs) {
    return first== rhs.first && second==rhs.second;
  }

  bool operator!=(const SimpleLor &rhs) {
    return    !((*this)==rhs);
  }
  struct Comparator { 
     int operator()(const SimpleLor &a, const SimpleLor &b) const {
      if(a.first<b.first) return 1;
      if(a.first>b.first) return 0;      
      return a.second<b.second;
    }
  };

  class Iterator ;

  static Iterator begin() ;
  static Iterator end() ;

  static const Comparator less;
private:
  static int n_lors_;
  static int n_detectors_;
};
   

class SimpleLor::Iterator {
 public:
 Iterator(int first, int second):lor(first,second){};

    SimpleLor::Iterator &operator++() {
      lor.second++;
      if(lor.second>=lor.first){
	lor.second=0;
	lor.first++;
      }
      return *this;
    }

    SimpleLor &operator*() {
      return lor;
    }

    bool operator!=(const SimpleLor::Iterator rhs) {
      return lor!=rhs.lor;
    }
 private:
    SimpleLor lor;
  };


inline SimpleLor::Iterator SimpleLor::begin() {return SimpleLor::Iterator(1,0);}
inline SimpleLor::Iterator SimpleLor::end() {return SimpleLor::Iterator(n_detectors_,0);}


std::ostream &operator<<(std::ostream &out,const SimpleLor &lor);
