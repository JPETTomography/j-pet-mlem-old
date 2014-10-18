#pragma once

namespace PET2D {
namespace Strip {
namespace GPU {

template <typename F> struct EventsSOA {
  F* z_u;
  F* z_d;
  F* dl;

  EventsSOA(const Event<F>* source, size_t n_events) {

    // temporary CPU side SOA
    F* cpu_z_u = new F[n_events];
    F* cpu_z_d = new F[n_events];
    F* cpu_dl = new F[n_events];

    // fill CPU side SOA
    for (int i = 0; i < n_events; ++i) {
      cpu_z_u[i] = source[i].z_u;
      cpu_z_d[i] = source[i].z_d;
      cpu_dl[i] = source[i].dl;
    }

    size_t mem_size = n_events * sizeof(F);

    // create CUDA side pointers
    cudaMalloc(&z_u, mem_size);
    cudaMalloc(&z_d, mem_size);
    cudaMalloc(&dl, mem_size);

    // copy CPU to device
    cudaMemcpy(z_u, cpu_z_u, mem_size, cudaMemcpyHostToDevice);
    cudaMemcpy(z_d, cpu_z_d, mem_size, cudaMemcpyHostToDevice);
    cudaMemcpy(dl, cpu_dl, mem_size, cudaMemcpyHostToDevice);

    delete[] cpu_z_u;
    delete[] cpu_z_d;
    delete[] cpu_dl;
  }

  ~EventsSOA() {
    cudaFree(z_u);
    cudaFree(z_d);
    cudaFree(dl);
  }
};
}  // GPU
}  // Strip
}  // PET2D
