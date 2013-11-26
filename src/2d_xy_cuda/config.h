#ifndef _CONFIG_H_
#define _CONFIG_H_

#define CUDA_CALLABLE_MEMBER __device__ __host__
#define THREAD_PER_PIXEL 1
#define THREAD_PER_EACH_PIXEL 0

#define NUMBER_OF_DETECTORS 64
#define LORS (NUMBER_OF_DETECTORS*(NUMBER_OF_DETECTORS + 1)) / 2

#define TP 6.28318530f
#define ITP 0.15915494f

#endif
