#pragma once

#include <iomanip>
#if _OPENMP
#include <omp.h>
#endif
#include <iostream>

namespace PET2D {
namespace Barrel {

/// Drives Monte-Carlo system matrix construction
template <typename DetectorType,
          typename MatrixType,
          typename FType = double,
          typename SType = int>
class MonteCarlo {
  using Detector = DetectorType;
  using Event = typename Detector::Event;
  using Matrix = MatrixType;
  using F = FType;
  using S = SType;
  using SS = typename std::make_signed<S>::type;
  using LOR = typename Matrix::LOR;

 public:
  MonteCarlo(const Detector& detector,
             Matrix& matrix,
             F pixel_size,
             S start_pixel = static_cast<S>(0))
      : detector(detector),
        matrix(matrix),
        pixel_size(pixel_size),
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

    bool tof = (tof_step > static_cast<F>(0));
    util::random::uniform_real_distribution<F> one_dis(0, 1);
    util::random::uniform_real_distribution<F> phi_dis(0, F(M_PI));

    matrix.add_emissions(n_emissions);

#if _OPENMP && !_MSC_VER
// We need to try catch inside OpenMP thread, otherwise we will not see the
// error thrown.
#define TRY try {
#define CATCH                     \
  }                               \
  catch (std::string & ex) {      \
    std::cerr << ex << std::endl; \
    throw(ex);                    \
  }
#else
#define TRY
#define CATCH
#endif

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
      TRY;
      auto pixel = matrix.pixel_at_index(i_pixel);

      if (pixel.x < start_pixel || pixel.y < start_pixel ||
          (pixel.x * pixel.x + pixel.y * pixel.y) * pixel_size * pixel_size >
              detector.barrel.fov_radius * detector.barrel.fov_radius)
        continue;
#if 0
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
        typename DetectorType::Response response;

        Event event(rx, ry, angle);
        auto hits = detector.detect(l_gen, model, event, response);

        // do we have hit on both sides?
        if (hits >= 2) {
          if (o_collect_mc_matrix) {
            if (response.lor.first == response.lor.second) {
              std::ostringstream msg;
              msg << __FUNCTION__ << " invalid LOR in Monte-Carlo ("
                  << response.lor.first << ", " << response.lor.second << ")";
              throw(msg.str());
            }
            matrix.hit_lor(response.lor, 0, i_pixel, 1);
          }

          if (o_collect_pixel_stats) {
            matrix.hit(i_pixel);
          }
          pixel_hit_count++;

        }  // if (hits>=2)
      }    // loop over emmisions from pixel
#endif
      matrix.compact_pixel_index(i_pixel);
      CATCH;
    }
  }

 private:
  const Detector& detector;
  Matrix& matrix;
  F pixel_size;
  F tof_step;
  S start_pixel;
};
}  // Barrel
}  // PET2D
