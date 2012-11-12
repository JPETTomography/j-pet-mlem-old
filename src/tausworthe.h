#ifndef _tausworthe_
#define _tausworthe_

typedef  double FLOAT ;

inline unsigned taus_step(unsigned &z, int S1, int S2, int S3, unsigned M) {
    unsigned b = (((z << S1)^z)>>S2);
    return z = (((z &M) << S3)^b);
  }

inline   unsigned LCG_step(unsigned &z, unsigned A, unsigned C) {
    return z = (A*z+C);
  }
inline   FLOAT rnd(unsigned *z) {
  return 2.328306436538696e-10*(FLOAT)(
          taus_step(z[0], 13, 19, 12, 4294967294u)^
          taus_step(z[1], 2, 25, 4, 4294967288u)^
          taus_step(z[2], 3, 11, 17, 4294967280u)^
          LCG_step(z[3], 1664525, 1013904223u) );
    }
inline  FLOAT rnd_nbc(unsigned *z1, unsigned *z2, unsigned *z3, unsigned *z4) {

  return 2.3283064e-10*(
      taus_step(*z1, 13, 19, 12, 4294967294u)^
      taus_step(*z2, 2, 25, 4, 4294967288u)^
      taus_step(*z3, 3, 11, 17, 4294967280u)^
      LCG_step(*z4, 1664525, 1013904223u) );
}

//    return 2.328306436538696e-10*(

inline   FLOAT rnd_np(unsigned  &z1, unsigned &z2, unsigned &z3, unsigned &z4) {

  return 2.328306436538696e-10*(FLOAT)(
          taus_step(z1, 13, 19, 12, 4294967294ul)^
          taus_step(z2, 2, 25, 4, 4294967288ul)^
          taus_step(z3, 3, 11, 17, 4294967280ul)^
          LCG_step(z4, 1664525, 1013904223ul) );
}
#endif
