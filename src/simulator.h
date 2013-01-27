#pragma once

#include "matrix_monte_carlo.h"
#include "detector_ring.h"

template <typename DetectorRingType, typename SystemMatrixType, typename F =
              double>
class Simulator {
  typedef typename DetectorRingType::LOR LOR;

 public:
  Simulator(DetectorRingType& detector_ring, SystemMatrixType& system_matrix)
      : detector_ring_(detector_ring),
        system_matrix_(system_matrix) {
  }

  /// Executes Monte-Carlo system matrix generation for given detector ring
  /// @param gen   random number generator
  /// @param model acceptance model (returns bool for call operator with given length)
  template <typename RandomGenerator, typename AcceptanceModel>
  void mc(RandomGenerator& gen,
          AcceptanceModel model,
          size_t n_mc_emissions,
          bool o_collect_mc_matrix = true,
          bool o_collect_pixel_stats = true) {

    uniform_real_distribution<> one_dis(0., 1.);
    uniform_real_distribution<> phi_dis(0., M_PI);

    int n_pixels_2 = system_matrix_.get_n_pixels() / 2;
    F s_pixel = system_matrix_.pixel_size();
    system_matrix_.increase_n_emissions(n_mc_emissions);

#if _OPENMP
    // OpenMP uses passed random generator as seed source for
    // thread local random generators
    RandomGenerator mp_gens[omp_get_max_threads()];
    for (auto t = 0; t < omp_get_max_threads(); ++t) {
      mp_gens[t].seed(gen());
    }

#pragma omp parallel for schedule(dynamic)
#endif
    // iterating only triangular matrix,
    // being upper right part or whole system matrix
    // descending, since biggest chunks start first, but may end last
    for (ssize_t y = n_pixels_2 - 1; y >= 0; --y) {
      for (auto x = 0; x <= y; ++x) {

        if ((x * x + y * y) * s_pixel * s_pixel >
            detector_ring_.fov_radius() * detector_ring_.fov_radius())
          continue;

        for (auto n = 0; n < n_mc_emissions; ++n) {
#if _OPENMP
          auto& l_gen = mp_gens[omp_get_thread_num()];
#else
          auto& l_gen = gen;
#endif
          auto rx = (x + one_dis(l_gen)) * s_pixel;
          auto ry = (y + one_dis(l_gen)) * s_pixel;

          // ensure we are within a triangle
          if (rx > ry)
            continue;

          auto angle = phi_dis(l_gen);
          LOR lor;
          auto hits =
              detector_ring_.emit_event(l_gen, model, rx, ry, angle, lor);

          // do we have hit on both sides?
          if (hits >= 2) {
            auto i_pixel = system_matrix_.t_pixel_index(x, y);

            if (o_collect_mc_matrix) {
              system_matrix_.add_to_t_matrix(lor, i_pixel);
            }

            if (o_collect_pixel_stats) {
              system_matrix_.add_hit(i_pixel);
            }
          }  //if (hits>=2)
        }    // loop over emmisions from pixel
      }
    }
  }

 private:
  DetectorRingType& detector_ring_;
  SystemMatrixType& system_matrix_;

  int n_emissions_;

};
