#ifndef SENSITIVITY_MAPPER
#define SENSITIVITY_MAPPER

#if _OPENMP
#include <omp.h>
#endif

#include "util/random.h"
#include "3d/geometry/voxel_set.h"
#include "3d/geometry/event_generator.h"

namespace PET3D {
namespace Hybrid {
template <typename Scanner> class SensitivityMapper {
 public:
  using F = typename Scanner::F;
  using S = typename Scanner::S;
  using VoxelSet = PET3D::VoxelSet<F, S>;
  using Event = typename Scanner::Event;

  SensitivityMapper(Scanner& scanner, VoxelSet& voxel_set)
      : scanner(scanner), voxel_set(voxel_set), one_dis(0, 1){};

  template <typename RNG, typename AcceptanceModel>
  void map(int i_voxel,
           const Voxel<S>& voxel,
           RNG& gen,
           AcceptanceModel& model,
           int n_emissions) {
    auto pixel_size = voxel_set.v_grid.pixel_grid.pixel_size;

    PET3D::Point<F> ll =
        voxel_set.v_grid.lower_left_at(voxel.ix, voxel.iy, voxel.iz);

    // std::cout<<"emitting from pixel at "<<ll.x<<" "<<ll.y<<" "<<ll.z<<"\n";

    for (int i = 0; i < n_emissions; ++i) {

      F rx = ll.x + one_dis(gen) * pixel_size;
      F ry = ll.y + one_dis(gen) * pixel_size;
      F rz = ll.z + one_dis(gen) * pixel_size;


      auto dir = direction(gen);
      //std::cout<<dir.x<<" "<<dir.y<<" "<<dir.z<<"\n";
      Event event(PET3D::Point<F>(rx, ry, rz), dir);

      typename Scanner::Response response;
      auto hits = scanner.detect(gen, model, event, response);
      if (hits >= 2) {
        voxel_set.value(i_voxel) += F(1.0);
      }
    }
  }

  template <typename RNG, typename AcceptanceModel>
  void map(RNG gen, AcceptanceModel& model, int n_emissions) {
#if _OPENMP
    // OpenMP uses passed random generator as seed source for
    // thread local random generators
    RNG* mp_gens = new (alloca(sizeof(RNG) * omp_get_max_threads()))
        RNG[omp_get_max_threads()];
    for (auto t = 0; t < omp_get_max_threads(); ++t) {
      mp_gens[t].seed(gen());
    }

#pragma omp parallel for schedule(dynamic)
// #pragma omp parallel for
#endif
    for (int i_voxel = 0; i_voxel < voxel_set.size(); ++i_voxel) {
#if _OPENMP
      auto& l_gen = mp_gens[omp_get_thread_num()];
#else
      auto& l_gen = gen;
#endif
      auto voxel = voxel_set.voxel(i_voxel);

      map(i_voxel, voxel, l_gen, model, n_emissions);
    }
  }

 private:
  Scanner& scanner;
  VoxelSet& voxel_set;
  util::random::uniform_real_distribution<F> one_dis;
  SphericalDistribution<F> direction;
};
}
}

#endif  // SENSITIVITY_MAPPER
