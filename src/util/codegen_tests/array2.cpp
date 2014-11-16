//  c++ -O2 -fomit-frame-pointer -std=c++11 -c -S -o - array2.cpp | c++filt
//
//  add(float, float):
//  	addss	%xmm1, %xmm0
//  	retq
//
//  nvcc -std=c++11 -x cu -ptx -o - array2.cpp | sed s,_param_,:, | c++filt -n
//
//  THIS IS EXPECTED, BUT AS OF CUDA 7.0 EA DOES NOT WORK!
//
//  .visible .func  (.param .b32 func_retval0) add(int, int)(
//  	.param .b32 add(int, int):0,
//  	.param .b32 add(int, int):1
//  )
//  {
//  	.reg .s32 	%r<4>;
//
//
//  	ld.param.u32 	%r1, [add(int, int):0];
//  	ld.param.u32 	%r2, [add(int, int):1];
//  	add.s32 	%r3, %r2, %r1;
//  	st.param.b32	[func_retval0+0], %r3;
//  	ret;
//  }

#define _ __device__
#if !__CUDACC__
#define __device__
#endif
#include "../array.h"

struct Bar {
  float foo, bar;
  _ Bar(float foo, float bar) : foo(foo), bar(bar) {}
};

_ float add(float a, float b) {
  util::array<4, Bar> bars;
  Bar bar(a, b);
  bars.push_back(bar);
  return bars[0].foo + bars[0].bar;
}
