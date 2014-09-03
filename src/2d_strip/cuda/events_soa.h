#pragma once

template <typename F> struct Events_SOA {
  F* z_u;
  F* z_d;
  F* dl;
};

template <typename F> Events_SOA<F> cuda_malloc_events_soa(size_t n_events) {
  Events_SOA<F> soa;
  size_t mem_size = n_events * sizeof(F);

  cudaMalloc(&soa.z_u, mem_size);
  cudaMalloc(&soa.z_d, mem_size);
  cudaMalloc(&soa.dl, mem_size);

  return soa;
}

template <typename F> Events_SOA<F> malloc_events_soa(size_t n_events) {
  Events_SOA<F> soa;
  soa.z_u = new F[n_events];
  soa.z_d = new F[n_events];
  soa.dl = new F[n_events];
  return soa;
}

template <typename F> void free_events_soa(Events_SOA<F> events) {
  delete[] events.z_u;
  delete[] events.z_d;
  delete[] events.dl;
}

template <typename F>
void transform_events_aos_to_soa(Events_SOA<F> dest,
                                 const Event<F>* source,
                                 size_t n_events) {

  for (int i = 0; i < n_events; ++i) {
    dest.z_u[i] = source[i].z_u;
    dest.z_d[i] = source[i].z_d;
    dest.dl[i] = source[i].dl;
  }
}

template <typename F>
void copy_events_soa_to_device(Events_SOA<F> dest,
                               Events_SOA<F> source,
                               size_t n_events) {
  size_t mem_size = n_events * sizeof(F);
  cudaMemcpy(dest.z_u, source.z_u, mem_size, cudaMemcpyHostToDevice);
  cudaMemcpy(dest.z_d, source.z_d, mem_size, cudaMemcpyHostToDevice);
  cudaMemcpy(dest.dl, source.dl, mem_size, cudaMemcpyHostToDevice);
}

template <typename F>
void load_events_to_gpu(const Event<F>* events,
                        Events_SOA<F> gpu_events,
                        size_t n_events) {
  Events_SOA<F> cpu_events = malloc_events_soa<F>(n_events);
  transform_events_aos_to_soa(cpu_events, events, n_events);
  copy_events_soa_to_device(gpu_events, cpu_events, n_events);
  free_events_soa(cpu_events);
}
