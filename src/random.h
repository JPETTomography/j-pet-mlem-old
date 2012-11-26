#pragma once

#include <cstdlib>

class tausworthe {
public:


  typedef unsigned int result_type;
  typedef long     int seed_type;

  static result_type min()  {return 0;}
  static result_type max()  {
    return std::numeric_limits<result_type>::max();
  }

  tausworthe(seed_type a_seed = 121245) {
    seed(a_seed);
  }

  result_type operator () () {
    taus_step(seeds[0], 13, 19, 12, 4294967294u);
    taus_step(seeds[1],  2, 25,  4, 4294967288u);
    taus_step(seeds[2],  3, 11, 17, 4294967280u);
    LCG_step (seeds[3],   1664525u, 1013904223u);
    return seeds[0]^seeds[1]^seeds[2]^seeds[3];
  }

  void seed(seed_type a_seed) {
    srand48(a_seed);
    for(int i = 0; i < 4; ++i) {
      result_type r;
      while( (r = lrand48()) < 128);
      seeds[i] = r;
    }
  }

private:
  unsigned int seeds[4];

  template<typename T>
  void taus_step(T &z, int S1, int S2, int S3, T M) {
    unsigned b = (((z << S1)^z)>>S2);
    z = (((z &M) << S3)^b);
  }

  template<typename T>
  void LCG_step(T &z, T A, T C) {
    z = (A*z+C);
  }
};

template<typename F = double>
class uniform_real_distribution {


public:
  typedef F result_type;

  uniform_real_distribution(F a = 0., F b = 1.)
  : s(b - a)
  , o(a)
  {}

  template<class Generator>
  result_type operator()(Generator& g) {
    return g() * size() * scale<Generator>() + offset()-
      static_cast<F>(Generator::min())/range<Generator>();
  }

private:
  result_type size()   const { return s; }
  result_type offset() const { return o; }



  template<class Generator>
  static result_type range() {
   return static_cast<F>(Generator::max())-static_cast<F>(Generator::min());
  }

  template<class Generator>
  static result_type scale() {
    return static_cast<F>(1.0)/range<Generator>();
  }

  F s, o;
};
