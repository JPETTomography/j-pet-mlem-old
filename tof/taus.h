#ifndef __TAUS_H__
#define __TAUS_H__

#include <cstdlib>

class taus_array {

public:

  taus_array(int n):n_generators_(n) {
    seeds_ = new unsigned[4*n_generators_];
    gen_seeds(121245);
  }

inline  double rnd(int i) {
  taus_step(seeds_[i*4], 13, 19, 12, 4294967294u);
  taus_step(seeds_[i*4+1], 2, 25, 4, 4294967288u);
  taus_step(seeds_[i*4+2], 3, 11, 17, 4294967280u);
  LCG_step(seeds_[i*4+3], 1664525u, 1013904223u);
  return 2.328306436538696e-10*(double)(seeds_[i*4]^
          seeds_[i*4+1]^
          seeds_[i*4+2]^
          seeds_[i*4+3]
          );
  }

inline  unsigned int  irnd(int i) {
  taus_step(seeds_[i*4], 13, 19, 12, 4294967294u);
  taus_step(seeds_[i*4+1], 2, 25, 4, 4294967288u);
  taus_step(seeds_[i*4+2], 3, 11, 17, 4294967280u);
  LCG_step(seeds_[i*4+3], 1664525u, 1013904223u);
  return (seeds_[i*4]^
    seeds_[i*4+1]^
    seeds_[i*4+2]^
    seeds_[i*4+3]
    );
  }

  void gen_seeds(long int);
  void set_seeds(const unsigned *seed);
  void get_seeds(unsigned *seed) const;

private:
  int n_generators_;
  unsigned *seeds_;

void taus_step(unsigned &z, int S1, int S2, int S3, unsigned M) {
    unsigned b = (((z << S1)^z)>>S2);
    z = (((z &M) << S3)^b);
  }

void LCG_step(unsigned &z, unsigned A, unsigned C) {
    z = (A*z+C);
  }
};

void taus_array::gen_seeds(long seed) {
  srand48(seed);
  for(int i = 0;i<4*n_generators_;++i) {
    unsigned r;
    while( (r = lrand48()) < 128) {};
    seeds_[i] = r;
  }
}

void taus_array::set_seeds(const unsigned   *seed) {
  for(int i = 0;i<4*n_generators_;++i) {
    seeds_[i] = seed[i];
  }
}

void taus_array::get_seeds(unsigned   *seed) const {
  for(int i = 0;i<4*n_generators_;++i) {
    seed[i] = seeds_[i];
  }
}
#endif
