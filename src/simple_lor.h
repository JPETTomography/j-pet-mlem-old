#pragma once

class SimpleLor {
public:
  int first;
  int second;

  static int t_index(const SimpleLor &lor) {
    SimpleLor ordered_lor=lor;
    if(ordered_lor.first< ordered_lor.second)
      swap(ordered_lor.first, ordered_lor.second);
    return (ordered_lor.first*(ordered_lor.first+1))/2+ordered_lor.second;
  }

  static void init(int n_detectors_a) {
    n_detectors_=n_detectors_a;
    n_lors_     =(n_detectors*(n_detectors-1))/2;
  }

  static const Comparator less;
private:
  static int n_lors_;
  static int n_detectors_;

  struct Comparator {
    int operator()(const SimpleLor &a, const SimpleLor &b) {
      if(a.first<b.first) return 1;
      if(a.first>b.first) return 0;      return a.second<b.second;
    }
  };
};
