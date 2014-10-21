#pragma once
#include <iomanip>
#if _OPENMP
#include <omp.h>
#endif

namespace PET2D {
namespace Barrel {

/// Drives Monte-Carlo system matrix construction
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
             F tof_step,
             S start_pixel = static_cast<S>(0))
      : detector_ring(detector_ring),
        matrix(matrix),
        pixel_size(pixel_size),
        tof_step(tof_step),
        start_pixel(start_pixel) {}

  /// Executes Monte-Carlo system matrix generation for given detector ring
  template <typename RandomGenerator,
            typename AcceptanceModel,
            typename ProgressCallback>
  void operator()(
      RandomGenerator& gen,              ///< random number generator
      AcceptanceModel model,             ///< acceptance model
      S n_emissions,                     ///< number of emissions generated
      ProgressCallback progress,         ///< progress callback
      bool o_collect_mc_matrix = true,   ///< enable matrix generation
      bool o_collect_pixel_stats = true  ///< enable pixel stats
      ) {
    if (n_emissions <= 0)
      return;

    auto n_positions = detector_ring.n_positions(tof_step, model.max_bias());
    bool tof = (tof_step > static_cast<F>(0));
    util::random::uniform_real_distribution<F> one_dis(0, 1);
    util::random::uniform_real_distribution<F> phi_dis(0, F(M_PI));

    matrix.add_emissions(n_emissions);

#if _OPENMP
    // OpenMP uses passed random generator as seed source for
    // thread local random generators
    RandomGenerator* mp_gens =
        new (alloca(sizeof(RandomGenerator) * omp_get_max_threads()))
            RandomGenerator[omp_get_max_threads()];
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
      progress(i_pixel);

      auto pixel = matrix.pixel_at_index(i_pixel);

      if (pixel.x < start_pixel || pixel.y < start_pixel ||
          (pixel.x * pixel.x + pixel.y * pixel.y) * pixel_size * pixel_size >
              detector_ring.fov_radius() * detector_ring.fov_radius())
        continue;

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
        F position = 0;
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
              msg << __FUNCTION__ << " invalid LOR in Monte-Carlo ("
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

#ifdef GPU_TOF_TEST

  template <typename RandomGenerator, typename AcceptanceModel>
  void test(RandomGenerator& gen, AcceptanceModel model, S n_emissions) {
    if (n_emissions <= 0)
      return;

    auto n_positions = detector_ring.n_positions(tof_step, model.max_bias());
    bool tof = (tof_step > static_cast<F>(0));
    uniform_real_distribution<F> one_dis(0, 1);
    uniform_real_distribution<F> phi_dis(0, F(M_PI));

    matrix.add_emissions(n_emissions);

    for (auto i_pixel = 0; i_pixel < matrix.total_n_pixels_in_triangle();
         ++i_pixel) {

      auto pixel = matrix.pixel_at_index(i_pixel);

      printf("PIXEL(%d,%d)\n", pixel.x, pixel.y);
      printf("N_POSITIONS: %d\n", n_positions);

      if ((pixel.x * pixel.x + pixel.y * pixel.y) * pixel_size * pixel_size >
          detector_ring.fov_radius() * detector_ring.fov_radius())
        continue;

      for (auto n = 0; n < n_emissions; ++n) {

        auto& l_gen = gen;

        auto rx = pixel.x * pixel_size;
        auto ry = pixel.y * pixel_size;

        for (auto angle = 0.0; angle < 3.14; angle += 0.1) {
          LOR lor;
          F position = 0;
          auto hits = detector_ring.emit_event(
              l_gen, model, rx, ry, angle, lor, position);

          S quantized_position = 0;
          if (tof)
            quantized_position = detector_ring.quantize_position(
                position, tof_step, n_positions);

          if (hits > 1) {
            printf("TOF: %d Position: %f \n", quantized_position, position);
          }
          // matrix.hit_lor(lor, quantized_position, i_pixel, 1);
        }
      }  // loop over emmisions from pixel
    }
  }

#endif

 private:
  DetectorRing& detector_ring;
  Matrix& matrix;
  F pixel_size;
  F tof_step;
  S start_pixel;
};
}  // Barrel
}  // PET2D
