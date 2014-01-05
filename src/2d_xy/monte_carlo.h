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
      : detector_ring(detector_ring),
        matrix(matrix),
        pixel_size(pixel_size),
        tof_step(tof_step) {}

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

    auto n_positions = detector_ring.n_positions(tof_step, model.max_bias());
    bool tof = (tof_step > static_cast<F>(0));
    uniform_real_distribution<> one_dis(0., 1.);
    uniform_real_distribution<> phi_dis(0., M_PI);

    matrix.add_emissions(n_emissions);

#if _OPENMP
    // OpenMP uses passed random generator as seed source for
    // thread local random generators
    RandomGenerator mp_gens[omp_get_max_threads()];
    for (auto t = 0; t < omp_get_max_threads(); ++t) {
      mp_gens[t].seed(gen());
    }

#pragma omp parallel for schedule(dynamic)
// #pragma omp parallel for
#endif
    // iterating only triangular matrix,
    // being upper right part or whole system matrix
    // NOTE: we must iterate pixel indices instead of x, y since we need proper
    // thread distribution when issuing on MIC
    for (auto i_pixel = 0; i_pixel < matrix.total_n_pixels_in_triangle();
         ++i_pixel) {

      auto pixel = matrix.pixel_at_index(i_pixel);

      if ((pixel.x * pixel.x + pixel.y * pixel.y) * pixel_size * pixel_size >
          detector_ring.fov_radius() * detector_ring.fov_radius())
        continue;
#if _OPENMP
        auto thread_id = omp_get_thread_num();
#else
        auto thread_id = 0;
#endif

      std::cerr<<"thread "<<thread_id<<" starting pixel "<<pixel.x<<" "<<pixel.y<<std::endl;

      int pixel_hit_count = 0;
      for (auto n = 0; n < n_emissions; ++n) {
#if _OPENMP
        auto& l_gen = mp_gens[omp_get_thread_num()];
#else
        auto& l_gen = gen;
#endif
        auto rx = (pixel.x + one_dis(l_gen)) * pixel_size;
        auto ry = (pixel.y + one_dis(l_gen)) * pixel_size;

        // ensure we are within a triangle
        if (rx > ry)
          continue;

        auto angle = phi_dis(l_gen);
        LOR lor;
        F position = (F)0.0;
        auto hits = detector_ring.emit_event(
            l_gen, model, rx, ry, angle, lor, position);

        S quantized_position = 0;
        if (tof)
          quantized_position =
              detector_ring.quantize_position(position, tof_step, n_positions);
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
            matrix.hit_lor(lor, quantized_position, i_pixel, 1);
          }

          if (o_collect_pixel_stats) {
            matrix.hit(i_pixel);
          }
          pixel_hit_count++;

        }  // if (hits>=2)
      }    // loop over emmisions from pixel
      matrix.compact_pixel_index(i_pixel);
    }
  }

 private:
  DetectorRing& detector_ring;
  Matrix& matrix;
  F pixel_size;
  F tof_step;
};
