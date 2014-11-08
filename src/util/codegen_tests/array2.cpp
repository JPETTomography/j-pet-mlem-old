// run: c++ -O3 -fomit-frame-pointer -std=c++11 -c array2.cpp
// run: otool -tv array2.o
// expected:
//   array2.o:
//   (__TEXT,__text) section
//   __Z3addff:
//   0000000000000000	addss	%xmm1, %xmm0
//   0000000000000004	retq

#include "../array.h"

struct Bar {
  float foo, bar;
  Bar(float foo, float bar) : foo(foo), bar(bar) {}
};

float add(float a, float b) {
  util::array<4, Bar> bars;
  Bar bar(a, b);
  bars.push_back(bar);
  return bars[0].foo + bars[0].bar;
}
