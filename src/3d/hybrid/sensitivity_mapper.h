#pragma once

#if _OPENMP
#include <omp.h>
#endif

#include "util/random.h"
#include "3d/geometry/voxel_set.h"
#include "3d/geometry/event_generator.h"

namespace PET3D {
namespace Hybrid {
template <class ScannerClass> class SensitivityMapper {
 public:
  using Scanner = ScannerClass;
  using F = typename Scanner::F;
  using S = typename Scanner::S;
  using VoxelSet = PET3D::VoxelSet<F, S>;
  using Event = typename Scanner::Event;

  SensitivityMapper(Scanner& scanner, VoxelSet& voxel_set)
      : scanner(scanner), voxel_set(voxel_set), one_dis(0, 1) {}

  template <class RNG, typename AcceptanceModel>
  void map(int i_voxel,
           const Voxel<S>& voxel,
           RNG& rng,
           AcceptanceModel& model,
           int n_emissions) {
    auto pixel_size = voxel_set.v_grid.pixel_grid.pixel_size;

    PET3D::Point<F> ll =
        voxel_set.v_grid.lower_left_at(voxel.ix, voxel.iy, voxel.iz);

#if DEBUG
    std::cout << "emitting from pixel at " << ll.x << " " << ll.y << " " << ll.z
              << std::endl;
#endif

    for (int i = 0; i < n_emissions; ++i) {

      F rx = ll.x + one_dis(rng) * pixel_size;
      F ry = ll.y + one_dis(rng) * pixel_size;
      F rz = ll.z + one_dis(rng) * pixel_size;

      auto dir = direction(rng);
#if DEBUG
      std::cout << dir.x << " " << dir.y << " " << dir.z << std::endl;
#endif
      Event event(PET3D::Point<F>(rx, ry, rz), dir);

      typename Scanner::Response response;
      auto hits = scanner.detect(rng, model, event, response);
      if (hits >= 2) {
        voxel_set.value(i_voxel) += F(1.0);
      }
    }
  }

  template <class RNG, typename AcceptanceModel>
  void map(RNG& rng, AcceptanceModel& model, int n_emissions) {
#if _OPENMP
    // OpenMP uses passed random generator as seed source for
    // thread local random generators
    RNG* mp_rngs = new (alloca(sizeof(RNG) * omp_get_max_threads()))
        RNG[omp_get_max_threads()];
    for (auto t = 0; t < omp_get_max_threads(); ++t) {
      mp_rngs[t].seed(rng());
    }

#pragma omp parallel for schedule(dynamic)
// #pragma omp parallel for
#endif
    for (int i_voxel = 0; i_voxel < (int)voxel_set.size(); ++i_voxel) {
#if _OPENMP
      auto& l_rng = mp_rngs[omp_get_thread_num()];
#else
      auto& l_rng = rng;
#endif
      auto voxel = voxel_set.voxel(i_voxel);

      map(i_voxel, voxel, l_rng, model, n_emissions);
    }
  }

 private:
  Scanner& scanner;
  VoxelSet& voxel_set;
  util::random::uniform_real_distribution<F> one_dis;
  SphericalDistribution<F> direction;
};

}  // Hybrid
}  // PET3D
