#pragma once

#include "matrix_lor_major.h"
#include "detector_ring.h"
#include "pixel.h"

template <typename DetectorRingType,
          typename SystemMatrixType,
          typename FType = double,
          typename SType = int>
class MonteCarlo {
  typedef DetectorRingType DetectorRing;
  typedef SystemMatrixType SystemMatrix;
  typedef FType F;
  typedef SType S;
  typedef typename std::make_signed<S>::type SS;
  typedef typename SystemMatrix::LOR LOR;

 public:
  MonteCarlo(
      DetectorRing& detector_ring, SystemMatrix& system_matrix, F pixel_size)
      : detector_ring_(detector_ring),
        system_matrix_(system_matrix),
        pixel_size_(pixel_size) {
  }

  /// Executes Monte-Carlo system matrix generation for given detector ring
  /// @param gen   random number generator
  /// @param model acceptance model (returns bool for call operator with given length)
  template <typename RandomGenerator, typename AcceptanceModel>
  void operator()(RandomGenerator& gen,
                  AcceptanceModel model,
                  S n_mc_emissions,
                  bool o_collect_mc_matrix = true,
                  bool o_collect_pixel_stats = true) {

    uniform_real_distribution<> one_dis(0., 1.);
    uniform_real_distribution<> phi_dis(0., M_PI);

    auto n_pixels_2 = system_matrix_.n_pixels_in_row() / 2;
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
    for (SS y = n_pixels_2 - 1; y >= 0; --y) {
      for (auto x = 0; x <= y; ++x) {

        if ((x * x + y * y) * pixel_size_ * pixel_size_ >
            detector_ring_.fov_radius() * detector_ring_.fov_radius())
          continue;

        auto i_pixel = Pixel<S>(x, y).index();

        for (auto n = 0; n < n_mc_emissions; ++n) {
#if _OPENMP
          auto& l_gen = mp_gens[omp_get_thread_num()];
#else
          auto& l_gen = gen;
#endif
          auto rx = (x + one_dis(l_gen)) * pixel_size_;
          auto ry = (y + one_dis(l_gen)) * pixel_size_;

          // ensure we are within a triangle
          if (rx > ry)
            continue;

          auto angle = phi_dis(l_gen);
          LOR lor;
          auto hits =
              detector_ring_.emit_event(l_gen, model, rx, ry, angle, lor);

          // do we have hit on both sides?
          if (hits >= 2) {
            if (o_collect_mc_matrix) {
              system_matrix_.add_to_t_matrix(lor, i_pixel);
            }

            if (o_collect_pixel_stats) {
              system_matrix_.add_hit(i_pixel);
            }
          }  // if (hits>=2)
        }    // loop over emmisions from pixel

        system_matrix_.compact_pixel_index(i_pixel);
      }
    }
  }

 private:
  DetectorRing& detector_ring_;
  SystemMatrix& system_matrix_;
  F pixel_size_;
};
