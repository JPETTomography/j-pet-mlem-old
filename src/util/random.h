#pragma once

#include <cstdlib>
#if _MSC_VER
#include <ctime>
#endif

class tausworthe {
 public:
  typedef unsigned int result_type;
  typedef long seed_type;

  static result_type min() { return 0; }
  static result_type max() { return std::numeric_limits<result_type>::max(); }

  tausworthe(seed_type a_seed = 121245) { seed(a_seed); }

  result_type operator()() {
    taus_step(seeds[0], 13, 19, 12, 4294967294u);
    taus_step(seeds[1], 2, 25, 4, 4294967288u);
    taus_step(seeds[2], 3, 11, 17, 4294967280u);
    LCG_step(seeds[3], 1664525u, 1013904223u);
    return seeds[0] ^ seeds[1] ^ seeds[2] ^ seeds[3];
  }

  void seed(seed_type a_seed) {
#if !_MSC_VER
    srand48(a_seed);
    for (int i = 0; i < 4; ++i) {
      result_type r;
      while ((r = static_cast<result_type>(lrand48())) < 128)
        ;
      seeds[i] = r;
    }
#else
    srand(time(NULL));
    for (int i = 0; i < 4; ++i) {
      result_type r;
      while ((r = static_cast<result_type>(rand())) < 128)
        ;
      seeds[i] = r;
    }
#endif
  }

 private:
  unsigned int seeds[4];

  template <typename T> void taus_step(T& z, int S1, int S2, int S3, T M) {
    unsigned b = (((z << S1) ^ z) >> S2);
    z = (((z & M) << S3) ^ b);
  }

  template <typename T> void LCG_step(T& z, T A, T C) { z = (A * z + C); }
};

template <typename FType = double> class uniform_real_distribution {
 public:
  typedef FType result_type;

  uniform_real_distribution(result_type a = 0, result_type b = 1)
      : size_(b - a), offset_(a) {}

  template <class Generator> result_type operator()(Generator& g) {
    return g() * size() * scale<Generator>() + offset() -
           static_cast<result_type>(Generator::min()) / range<Generator>();
  }

  result_type size() const { return size_; }
  result_type offset() const { return offset_; }

  template <class Generator> static result_type range() {
    return static_cast<result_type>(Generator::max()) -
           static_cast<result_type>(Generator::min());
  }

  template <class Generator> static result_type scale() {
    return static_cast<result_type>(1) / range<Generator>();
  }

 private:
  result_type size_, offset_;
};
