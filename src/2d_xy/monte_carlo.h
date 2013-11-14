#pragma once
#include <iomanip>
#if _OPENMP
#include <omp.h>
#endif

template <typename DetectorRingType,
          typename MatrixType,
          typename FType = double,
          typename SType = int>
class MonteCarlo {
  typedef DetectorRingType DetectorRing;
  typedef MatrixType Matrix;
  typedef FType F;
  typedef SType S;
  typedef typename std::make_signed<S>::type SS;
  typedef typename Matrix::LOR LOR;

 public:
  MonteCarlo(DetectorRing& detector_ring,
             Matrix& matrix,
             F pixel_size,
             F tof_step)
      : detector_ring_(detector_ring),
        matrix_(matrix),
        pixel_size_(pixel_size),
        tof_step_(tof_step) {}

  /// Executes Monte-Carlo system matrix generation for given detector ring
  /// @param gen   random number generator
  /// @param model acceptance model (returns bool for call operator with given
  /// length)
  template <typename RandomGenerator, typename AcceptanceModel>
  void operator()(RandomGenerator& gen,
                  AcceptanceModel model,
                  S n_emissions,
                  bool o_collect_mc_matrix = true,
                  bool o_collect_pixel_stats = true) {
    if (n_emissions <= 0)
      return;

    auto n_positions = detector_ring_.n_positions(tof_step_, model.max_bias());
    bool tof = (tof_step_ > static_cast<F>(0));
    uniform_real_distribution<> one_dis(0., 1.);
    uniform_real_distribution<> phi_dis(0., M_PI);

    auto n_pixels_2 = matrix_.n_pixels_in_row() / 2;
    matrix_.add_emissions(n_emissions);

#if _OPENMP
    // OpenMP uses passed random generator as seed source for
    // thread local random generators
    RandomGenerator mp_gens[omp_get_max_threads()];
    for (auto t = 0; t < omp_get_max_threads(); ++t) {
      mp_gens[t].seed(gen());
    }

#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
    // iterating only triangular matrix,
    // being upper right part or whole system matrix
    // descending, since biggest chunks start first, but may end last
    for (SS y = n_pixels_2 - 1; y >= 0; --y) {
      for (auto x = 0; x <= y; ++x) {
        int thread_id = 0;
#if _OPENMP
        thread_id = omp_get_thread_num();
#endif

        if ((x * x + y * y) * pixel_size_ * pixel_size_ >
            detector_ring_.fov_radius() * detector_ring_.fov_radius())
          continue;

        auto i_pixel = Pixel<S>(x, y).index();
        int pixel_hit_count = 0;
        for (auto n = 0; n < n_emissions; ++n) {
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
          F position = (F)0.0;
          auto hits = detector_ring_.emit_event(
              l_gen, model, rx, ry, angle, lor, position);

          S quantized_position = 0;
          if (tof)
            quantized_position = detector_ring_.quantize_position(
                position, tof_step_, n_positions);
#ifdef DEBUG
          std::cerr << "quantized_position " << quantized_position << std::endl;
#endif
          // do we have hit on both sides?
          if (hits >= 2) {
            if (o_collect_mc_matrix) {
              if (lor.first == lor.second) {
                std::ostringstream msg;
                msg << __PRETTY_FUNCTION__ << " invalid LOR in Monte-Carlo ("
                    << lor.first << ", " << lor.second << ")";
                throw(msg.str());
              }
              matrix_.hit_lor(lor, quantized_position, i_pixel, 1);
            }

            if (o_collect_pixel_stats) {
              matrix_.hit(i_pixel);
            }
            pixel_hit_count++;

          }  // if (hits>=2)
        }    // loop over emmisions from pixel
        matrix_.compact_pixel_index(i_pixel);
      }
    }
  }

 private:
  DetectorRing& detector_ring_;
  Matrix& matrix_;
  F pixel_size_;
  F tof_step_;
};
