// run: c++ -O3 -fomit-frame-pointer -std=c++11 -c array1.cpp
// run: otool -tv array1.o
// expected:
//   array1.o:
//   (__TEXT,__text) section
//   __Z3addii:
//   0000000000000000	addl	%esi, %edi
//   0000000000000002	movl	%edi, %eax
//   0000000000000004	retq

#include "../array.h"

struct Foo {
  int foo;
  Foo(int foo) : foo(foo) {}
};

int add(int a, int b) {
  util::array<2, Foo> foos{ a, b };
  return foos[0].foo + foos[1].foo;
}
