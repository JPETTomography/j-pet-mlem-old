//  c++ -O2 -fomit-frame-pointer -std=c++11 -c -S -o - array1.cpp | c++filt
//
//  add(int, int):
//  	addl	%esi, %edi
//  	movl	%edi, %eax
//  	retq
//
//  nvcc -std=c++11 -x cu -ptx -o - array1.cpp | sed s,_param_,:, | c++filt -n
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

struct Foo {
  int foo;
  _ Foo(int foo) : foo(foo) {}
};

_ int add(int a, int b) {
  util::array<4, Foo> foos{ a, b };
  return foos[0].foo + foos[1].foo;
}
