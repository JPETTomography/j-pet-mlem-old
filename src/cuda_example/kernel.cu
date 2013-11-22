// if we don't include that Qt Creator will show many errors
#include <cuda_runtime.h>
#include <sys/time.h>
#include <stdio.h>
#include "config.h"

/*-----------------------------------------GENERATORS-----------------------------------------------------*/

__constant__ unsigned int shift1[4] = { 6, 2, 13, 3 };
__constant__ unsigned int shift2[4] = { 13, 27, 21, 12 };
__constant__ unsigned int shift3[4] = { 18, 2, 7, 13 };
__constant__ unsigned int offset[4] = { 4294967294, 4294967288,
                                        4294967280, 4294967168 };

__shared__ unsigned int randStates[32];

//-----------end state rng
//---------------------------------------------------------//

__device__ float devData;

__device__ inline unsigned TausStep(unsigned& z,
                                    int S1,
                                    int S2,
                                    int S3,
                                    unsigned M) {
  unsigned b = (((z << S1) ^ z) >> S2);
  return z = (((z & M) << S3) ^ b);
}

__device__ inline unsigned LCGStep(unsigned& z, unsigned A, unsigned C) {
  return z = (A * z + C);
}

__device__ inline float HybridTaus(unsigned& z1,
                                   unsigned& z2,
                                   unsigned& z3,
                                   unsigned& z4) {
  return 2.3283064365387e-10f * (TausStep(z1, 14, 16, 15, 4294967294UL) ^
                                 TausStep(z2, 2, 44, 7, 4294967288UL) ^
                                 TausStep(z3, 3, 144, 17, 4294967280UL) ^
                                 LCGStep(z4, 1664525, 1013904223UL));
}

__device__ float rand_MWC_co(unsigned long long* x, unsigned int* a) {
  // Generate a random number [0,1)
  *x = (*x & 0xffffffffull) * (*a) + (*x >> 32);
  return __fdividef(
      __uint2float_rz((unsigned int)(*x)),
      (float)0x100000000);  // The typecast will truncate the x so that it is
                            // 0<=x<(2^32-1),__uint2float_rz ensures a round
                            // towards zero since 32-bit floating point cannot
                            // represent all integers that large. Dividing by
                            // 2^32 will hence yield [0,1)

}  // end __device__ rand_MWC_co

__device__ float rand_MWC_oc(unsigned long long* x, unsigned int* a) {
  // Generate a random number (0,1]
  return 1.0f - rand_MWC_co(x, a);
}  // end __device__ rand_MWC_oc

/*-------------------------------DATA STRUCTURES
 * DEFINITON-------------------------------------------------*/

struct data {

  float x, y;
  float f_data[8];
};

// x,y,z,w
//__builtin_align__(32)

struct Points {

  float x, y;
};

struct Hits {

  Points p[2];
};

struct Detectors {

  Points points[4];
};

struct Detector_Ring {

  Detectors detector_list[NUMBER_OF_DETECTORS];
};

struct Matrix_Element {

  float lor[LORS];
};

struct Secant_Points {

  float x1, y1, x2, y2;
};

struct Secant_Angle {

  float angle1, angle2;
};

struct Secant_Sections {

  int ss1, ss2;
};

/*****************************METHODS*************************************************/

CUDA_CALLABLE_MEMBER Secant_Points
secant(float x, float y, float angle, float radius) {

  float a = std::sin(angle);
  float b = -std::cos(angle);
  float c = a * x + b * y;

  // std::cout << "a: " << a << " b: " << b << " c: " << c << std::endl;

  // helper variables
  float b2 = b * b;
  float b2c = b2 * c;
  float ac = a * c;
  float a2_b2 = a * a + b2;
  float b_a2_b2 = b * a2_b2;

  float sq = sqrt(b2 * (-(c * c) + a2_b2 * radius * radius));
  float asq = a * sq;

  // std::cout << "sq: " << sq << " asq: " << asq << std::endl;

  Secant_Points secant_positions;

  secant_positions.x1 = (ac - sq) / a2_b2;
  secant_positions.y1 = ((b2c + asq) / b_a2_b2);
  secant_positions.x2 = (ac + sq) / a2_b2;
  secant_positions.y2 = ((b2c - asq) / b_a2_b2);

  return secant_positions;
}

CUDA_CALLABLE_MEMBER Secant_Angle secant_angles(Secant_Points& e) {

  Secant_Angle temp;
  temp.angle1 = atan2(e.y1, e.x1);
  temp.angle2 = atan2(e.y2, e.x2);

  return temp;
}

CUDA_CALLABLE_MEMBER int section(float angle, int n_detectors) {
  // converting angles to [0,2 Pi) interval
  float normalised_angle = angle > 0 ? angle : (float)2.0 * M_PI + angle;
  return static_cast<int>(round(normalised_angle * n_detectors * INV_TWO_PI)) %
         (n_detectors);
}

CUDA_CALLABLE_MEMBER Secant_Sections
secant_sections(Secant_Points& e, int n_detectors) {

  Secant_Angle angles = secant_angles(e);

  Secant_Sections temp;

  temp.ss1 = section(angles.angle1, n_detectors);
  temp.ss2 = section(angles.angle2, n_detectors);

  return temp;
}

CUDA_CALLABLE_MEMBER int intersections(float x,
                                       float y,
                                       float angle,
                                       Detector_Ring& ring,
                                       int detector_id,
                                       Hits& hit) {

  // float temp_x =
  // ring.detector_list[detector_id].get_element(3).get_locationx();
  // float temp_y =
  // ring.detector_list[detector_id].get_element(3).get_locationy();

  float p1_x = ring.detector_list[detector_id].points[3].x;
  float p1_y = ring.detector_list[detector_id].points[3].y;

  // std::cout << "temp_x: " <<  temp_x << " " << "temp_y: " << temp_y <<
  // std::endl;

  float a = std::sin(angle);
  float b = -std::cos(angle);
  float c = a * x + b * y;

  // std::cout << "a: " << a << " " << "b: " << b << " " << "c: " << c <<
  // std::endl;

  float v1 = a * p1_x + b * p1_y - c;

  // std::cout <<"V1: " <<  v1 << std::endl;

  int r = 0;

  for (int i = 0; i < 4; i++) {

    float p2_x = ring.detector_list[detector_id].points[i].x;
    float p2_y = ring.detector_list[detector_id].points[i].y;

    float v2 = a * p2_x + b * p2_y - c;

    if (v2 == 0.0f) {
      hit.p[r].x = ring.detector_list[detector_id].points[i].x;
      hit.p[r].y = ring.detector_list[detector_id].points[i].y;
      // v2 is crossing point
      r++;
      // std::cout << " " << p2_x << "  " << p2_y;
      if (r == 2)
        return r;
    } else if (v1 * v2 < 0.0f) {
      // calculate intersection

      float m = a * (p1_x - p2_x) + b * (p1_y - p2_y);
      /*
       std::cout
       << (c * (p1_x - p2_x) + b * (p2_x * p1_y - p1_x * p2_y)) / m
       << "  "
       << ((c * (p1_y - p2_y) + a * (p1_x * p2_y - p2_x * p1_y))
       / m) << std::endl;*/
      hit.p[r].x = (c * (p1_x - p2_x) + b * (p2_x * p1_y - p1_x * p2_y)) / m;
      hit.p[r].y = (c * (p1_y - p2_y) + a * (p1_x * p2_y - p2_x * p1_y)) / m;

      r++;

      if (r == 2)
        return r;
    }
    v1 = v2;
    p1_x = p2_x;
    p1_y = p2_y;
  }
  return r;
}

CUDA_CALLABLE_MEMBER float secant_angle(float x1, float y1) {
  return atan2(y1, x1);
}

CUDA_CALLABLE_MEMBER bool check_for_hits(int inner,
                                         int outer,
                                         int x,
                                         int y,
                                         float angle,
                                         int n_detectors,
                                         Detector_Ring& ring,
                                         int& detector,
                                         Hits& hit) {

  int points;

  int step = ((n_detectors + inner - outer) % n_detectors >
              (n_detectors + outer - inner) % n_detectors)
                 ? 1
                 : n_detectors - 1;
  int end = (outer + step) % n_detectors;
  for (int i = inner; i != end; i = (i + step) % n_detectors) {
    points = intersections(x, y, angle, ring, i, hit);

    if (points == 2) {

      detector = i;
      return true;
    }

    // check if we got 2 point intersection
    // then test the model against these points distance
    // if (points.size() == 2) {
    //   auto deposition_depth = model.deposition_depth(gen);
    //   if (deposition_depth < (points[1] - points[0]).length()) {
    //    detector = i;
    //    depth = deposition_depth;
    //    return true;
    //  }
  }

  return false;
}

/***************************************************************************************/

static cudaError err;

#define cuda(f, ...)                                        \
  if ((err = cuda##f(__VA_ARGS__)) != cudaSuccess) {        \
    fprintf(stderr, #f "() %s\n", cudaGetErrorString(err)); \
    exit(-1);                                               \
  }
#define cudathread_per_blockoSync(...) cuda(__VA_ARGS__)

double getwtime() {
  struct timeval tv;
  static time_t sec = 0;
  gettimeofday(&tv, NULL);
  if (!sec)
    sec = tv.tv_sec;
  return (double)(tv.tv_sec - sec) + (double)tv.tv_usec / 1e6;
}

__global__ void hello(char* a, int* b) {
  // increment chars with ints
  a[threadIdx.x] += b[threadIdx.x];
}

__device__ int lor_iterator(int& id1, int& id2) {

  if (id1 < id2) {
    int temp;
    temp = id2;
    id2 = id1;
    id1 = temp;
  }

  return ((id1 * (id1 + 1)) / 2) + id2;
}

__global__ void gpu_phantom_generation(int x,
                                       int y,
                                       int iteration,
                                       unsigned int* gpu_prng_seed,
                                       Matrix_Element* pixel_data,
                                       int threads,
                                       int pixels_in_row,
                                       float radius,
                                       float h_detector,
                                       float w_detector,
                                       float pixel_size) {

  int tid = ((blockIdx.x * blockDim.x) + threadIdx.x);

  unsigned int seed[4];

  seed[0] = gpu_prng_seed[4 * tid];
  seed[1] = gpu_prng_seed[4 * tid + 1];
  seed[2] = gpu_prng_seed[4 * tid + 2];
  seed[3] = gpu_prng_seed[4 * tid + 3];

  float inv_unit_prob_ = 0.1f;

  __shared__ Detector_Ring test_ring;

  Hits hit1;
  Hits hit2;

  float fov_radius = radius / M_SQRT2;

  if (threadIdx.x < NUMBER_OF_DETECTORS) {

    Detectors detector_base;

    detector_base.points[0].x =
        (w_detector / 2.0f) + radius + (h_detector / 2.0f);
    detector_base.points[0].y = h_detector / 2.0f;
    detector_base.points[1].x =
        (w_detector / 2.0f) + radius + (h_detector / 2.0f);
    detector_base.points[1].y = -h_detector / 2.0f;
    detector_base.points[2].x =
        (-w_detector / 2.0f) + radius + (h_detector / 2.0f);
    detector_base.points[2].y = -h_detector / 2.0f;
    detector_base.points[3].x =
        (-w_detector / 2.0) + radius + (h_detector / 2.0f);
    detector_base.points[3].y = h_detector / 2.0f;

    test_ring.detector_list[threadIdx.x] = detector_base;

    float angle = 2.0f * M_PI * threadIdx.x / NUMBER_OF_DETECTORS;
    float sin_phi = __sinf(angle);
    float cos_phi = __cosf(angle);

    for (int j = 0; j < 4; ++j) {

      float temp_x = test_ring.detector_list[threadIdx.x].points[j].x;
      float temp_y = test_ring.detector_list[threadIdx.x].points[j].y;

      test_ring.detector_list[threadIdx.x].points[j].x =
          temp_x * cos_phi - temp_y * sin_phi;
      test_ring.detector_list[threadIdx.x].points[j].y =
          temp_x * sin_phi + temp_y * cos_phi;
    }
  }

  __syncthreads();

  int detector1;
  int detector2;

  int i_inner;
  int i_outer;

#pragma unroll
  for (int i = 0; i < iteration; ++i) {

    if ((x * x + y * y) * pixel_size * pixel_size > fov_radius * fov_radius) {
      continue;
    }

    float rx =
        (x + HybridTaus(seed[0], seed[1], seed[2], seed[3])) * pixel_size;
    float ry =
        (y + HybridTaus(seed[0], seed[1], seed[2], seed[3])) * pixel_size;

    float angle = HybridTaus(seed[0], seed[1], seed[2], seed[3]) * M_PI;

    if (rx > ry) {
      continue;
    }

    // innetr and outer secant for circles
    Secant_Points inner_secant = secant(rx, ry, angle, radius);
    Secant_Points outer_secant = secant(rx, ry, angle, radius + h_detector);

    // hits per detector(if hits = 2 we got pair of detector, else generate
    // new random position and angle)

    i_inner = section(secant_angle(inner_secant.x1, inner_secant.y1),
                      NUMBER_OF_DETECTORS);
    i_outer = section(secant_angle(outer_secant.x1, inner_secant.y1),
                      NUMBER_OF_DETECTORS);

    if (!check_for_hits(i_inner,
                        i_outer,
                        x,
                        y,
                        angle,
                        NUMBER_OF_DETECTORS,
                        test_ring,
                        detector1,
                        hit1)) {
      continue;
    }

    i_inner = section(secant_angle(inner_secant.x2, inner_secant.y2),
                      NUMBER_OF_DETECTORS);
    i_outer = section(secant_angle(outer_secant.x2, inner_secant.y2),
                      NUMBER_OF_DETECTORS);

    if (!check_for_hits(i_inner,
                        i_outer,
                        x,
                        y,
                        angle,
                        NUMBER_OF_DETECTORS,
                        test_ring,
                        detector2,
                        hit2)) {
      continue;
    }

    float deposition_depth =
        -log(HybridTaus(seed[0], seed[1], seed[2], seed[3])) * inv_unit_prob_;

    // if(threadIdx.x == 1){
    // printf("TID:%d %d %d %d %d %f %f %f %f %f %f %f %f\n ", threadIdx.x, x,
    // y,
    //        detector1, detector2, hit1.p[0].x, hit1.p[0].y, hit1.p[1].x,
    //        hit1.p[1].y, hit2.p[0].x, hit2.p[0].y, hit2.p[1].x, hit2.p[1].y);
    // }
    // if(deposition_depth < (sqrt( (hit[1].x - hit[0].x) * (hit[1].x -
    // hit[0].x) + (hit[1].y - hit[0].y) * (hit[1].x - hit[0].x) ))) {

    atomicAdd(&pixel_data[blockIdx.x].lor[lor_iterator(detector1, detector2)],
              1.0f / iteration / threads);
    //}
  }

  gpu_prng_seed[4 * tid] = seed[0];
  gpu_prng_seed[4 * tid + 1] = seed[1];
  gpu_prng_seed[4 * tid + 2] = seed[2];
  gpu_prng_seed[4 * tid + 3] = seed[3];
}

void run_kernel(char* str, int* val, int str_size, int val_size) {

  // CUDA side variables
  char* cuda_str;
  int* cuda_val;

  cudaMalloc((void**)&cuda_str, str_size);
  cudaMalloc((void**)&cuda_val, val_size);
  cudaMemcpy(cuda_str, str, str_size, cudaMemcpyHostToDevice);
  cudaMemcpy(cuda_val, val, val_size, cudaMemcpyHostToDevice);

  dim3 dimBlock(10, 1);
  dim3 dimGrid(1, 1);
  // Qt Creator does not like that, but ask NVIDIA about fancy notation
  // hello << <dimGrid, dimBlock>>> (cuda_str, cuda_val);

  cudaMemcpy(str, cuda_str, str_size, cudaMemcpyDeviceToHost);

  cudaFree(cuda_str);
  cudaFree(cuda_str);
}

void mem_clean_lors(Matrix_Element* cpu_matrix, int number_of_blocks) {

  for (int i = 0; i < number_of_blocks; ++i) {
    for (int j = 0; j < LORS; ++j) {

      cpu_matrix[i].lor[j] = 0.f;
    }
  }
}

void phantom_kernel(int number_of_threads_per_block,
                    int number_of_blocks,
                    int n_emissions,
                    int pixels_in_row,
                    float radius,
                    float h_detector,
                    float w_detector,
                    float pixel_size) {

  dim3 blocks(number_of_blocks);
  dim3 threads(number_of_threads_per_block);

  unsigned int* cpu_prng_seed;

  cpu_prng_seed =
      (unsigned int*)malloc(number_of_blocks * number_of_threads_per_block * 4 *
                            sizeof(unsigned int));

  for (int i = 0; i < 4 * number_of_blocks * number_of_threads_per_block; ++i) {

    cpu_prng_seed[i] = 53445 + i;
  }

  int triangular_matrix_size =
      ((pixels_in_row / 2) * ((pixels_in_row / 2) + 1) / 2);

  Matrix_Element* cpu_matrix =
      (Matrix_Element*)malloc(number_of_blocks * sizeof(Matrix_Element));

  // unsigned int matrix_size = triangular_matrix_size * number_of_blocks;

  unsigned int* gpu_prng_seed;
  Matrix_Element* gpu_matrix_element;

  cuda(Malloc,
       (void**)&gpu_prng_seed,
       number_of_blocks * number_of_threads_per_block * 4 *
           sizeof(unsigned int));
  cuda(Malloc,
       (void**)&gpu_matrix_element,
       number_of_blocks * sizeof(Matrix_Element));

  cuda(
      Memcpy,
      gpu_prng_seed,
      cpu_prng_seed,
      number_of_blocks * number_of_threads_per_block * 4 * sizeof(unsigned int),
      cudaMemcpyHostToDevice);

  printf("GPU kernel start\n");
  printf(
      "Number of Detectors %d Numer of LORS: %d\n", NUMBER_OF_DETECTORS, LORS);

  double timer = getwtime();

  for (int j = pixels_in_row / 2 - 1; j >= 0; --j) {
    for (int i = 0; i <= j; ++i) {

      mem_clean_lors(cpu_matrix, number_of_blocks);

      cuda(Memcpy,
           gpu_matrix_element,
           cpu_matrix,
           number_of_blocks * sizeof(Matrix_Element),
           cudaMemcpyHostToDevice);

      printf("Pixel(%d,%d) n_emissions: %d \n", i, j, n_emissions);

      gpu_phantom_generation << <blocks, threads>>>
          (i,
           j,
           n_emissions,
           gpu_prng_seed,
           gpu_matrix_element,
           number_of_threads_per_block,
           pixels_in_row,
           radius,
           h_detector,
           w_detector,
           pixel_size);

      cudaThreadSynchronize();

      cuda(Memcpy,
           cpu_matrix,
           gpu_matrix_element,
           number_of_blocks * sizeof(Matrix_Element),
           cudaMemcpyDeviceToHost);

      for (int i = 0; i < LORS; i++) {
        float temp = 0.f;
        for (int j = 0; j < number_of_blocks; ++j) {

          temp += cpu_matrix[j].lor[i];
        }

        if (temp > 0.0f) {
          printf("%f\n", temp / number_of_blocks);
        }
      }
    }
  }
  double time = 0.0f;

  time = getwtime() - time;

  printf("time[s]: %f\n ", time);
  printf("time per pixel: %f\n", time / triangular_matrix_size);

  cuda(Free, gpu_prng_seed);
  cuda(Free, gpu_matrix_element);
}
