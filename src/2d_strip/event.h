#pragma once

#include <cmath>

#include "util/cuda/compat.h"
#include "util/utils.h"

template <typename F> struct Event {
  F z_u;
  F z_d;
  F dl;

  _ Event(F z_u, F z_d, F dl) : z_u(z_u), z_d(z_d), dl(dl) {}

  _ void transform(F R, F& tan, F& y, F& z) const {
    tan = this->tan(R);
    y = this->y(tan);
    z = this->z(y, tan);
  }

 private:
  _ F tan(const F R) const { return (z_u - z_d) / (2 * R); }

  _ F y(const F tan) const {
    return -F(0.5) * (dl / compat::sqrt(1 + tan * tan));
  }

  _ F z(const F y, const F tan) const {
    return F(0.5) * (z_u + z_d + (2 * y * tan));
  }
};

template <typename F> struct ImageSpaceEventTan;

template <typename F> struct ImageSpaceEventAngle {
  const F y;
  const F z;
  const F angle;

  _ ImageSpaceEventAngle(F y, F z, F angle) : y(y), z(z), angle(angle) {}

  _ ImageSpaceEventTan<F> to_tan() const {
    return ImageSpaceEventTan<F>(y, z, compat::tan(angle));
  }
};

template <typename F> struct ImageSpaceEventTan {
  const F y;
  const F z;
  const F tan;

  ImageSpaceEventTan(F y, F z, F tan) : y(y), z(z), tan(tan) {}

  _ ImageSpaceEventAngle<F> to_angle() const {
    return ImageSpaceEventAngle<F>(y, z, compat::atan(tan));
  }
};

template <typename FType> struct EllipseParameters {
  typedef FType F;
  F x, y, a, b;
  F angle;
  F n_emissions;

  _ EllipseParameters(F x, F y, F a, F b, F angle, F n_emissions)
      : x(x), y(y), a(a), b(b), angle(angle), n_emissions(n_emissions) {}
};

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
  size_t mem_size = n_events * sizeof(F);
  soa.z_u = (F*)safe_malloc(mem_size);
  soa.z_d = (F*)safe_malloc(mem_size);
  soa.dl = (F*)safe_malloc(mem_size);
  return soa;
}

template <typename F> void free_events_soa(Events_SOA<F> events) {
  free(events.z_u);
  free(events.z_d);
  free(events.dl);
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
