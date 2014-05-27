#include <2d_strip/event.h>

#pragma once

namespace gpu_config {

struct GPU_parameters {

  float R_distance, Scentilator_length, pixel_size;
  int n_pixels;
  float sigma;
  float dl;
  int number_of_blocks;
  int number_of_threads_per_block;
  int number_of_events;
  float inv_pow_sigma_dl;
  float inv_pow_sigma_z;
  float grid_size_y_;
  float grid_size_z_;
};
}

#define SOA_SIZE 180000000

template <typename T> struct soa_event {

  T z_u[SOA_SIZE];
  T z_d[SOA_SIZE];
  T dl[SOA_SIZE];

  void set_data_chunk(event<float>* data_chunk, int offset, int aos_data_size) {

    for (int i = 0; i < aos_data_size; i += offset) {

      for (int j = 0; j < offset; ++j) {

        if ((i + j) < aos_data_size) {

          z_u[i + j] = data_chunk[j].z_u;
          z_d[i + j] = data_chunk[j].z_d;
          dl[i + j] = data_chunk[j].dl;
        }
      }
    }
  }

  void set_data(event<float>* aos_data, int aos_data_size) {

    // data_size = aos_data_size;

    for (int i = 0; i < aos_data_size; ++i) {

      z_u[i] = aos_data[i].z_u;
      z_d[i] = aos_data[i].z_d;
      dl[i] = aos_data[i].dl;
    }
  }

  //    T* z_u;
  //    T* z_d;
  //    T* dl;
  //    int malloc_size;
  //    int data_size;
  /*
      soa_event(){

          printf("Allocate the SOA representation for event_list\n");
          z_u = (T*)malloc(N * sizeof(float));
          z_d = (T*)malloc(N * sizeof(float));
          dl = (T*)malloc(N * sizeof(float));
          malloc_size = N;
      }

      void set_data(event<float>* aos_data,int aos_data_size){

          data_size = aos_data_size;

          for(int i = 0; i < aos_data_size;++i){

              z_u[i] = aos_data[i].z_u;
              z_d[i] = aos_data[i].z_d;
              dl[i] = aos_data[i].dl;

          }
      }

      ~soa_event(){
          free(z_u);
          free(z_d);
          free(dl);
          printf("Free memory of SOA representation for event_list\n");
      }
  */
};
