#pragma once

#include <iomanip>
#if _OPENMP
#include <omp.h>
#endif
#include <iostream>

#include "3d/geometry/event_generator.h"
#include "3d/geometry/event.h"
#include "3d/geometry/point.h"

namespace PET3D {
namespace Hybrid {

/// Drives Monte-Carlo system matrix construction
template <class ScannerClass, class MatrixClass> class MonteCarlo {
  using Scanner = ScannerClass;
  using Event = typename Scanner::Event;
  using Matrix = MatrixClass;
  using F = typename Scanner::F;
  using S = typename Scanner::S;
  using SS = typename std::make_signed<S>::type;
  using LOR = typename Matrix::LOR;
  using Point = PET3D::Point<F>;

  static_assert(std::is_same<typename Matrix::S, S>::value,
                "matrix SType must be the same as detector SType");

 public:
  MonteCarlo(const Scanner& scanner,
             Matrix& matrix,
             F pixel_size,
             S start_pixel = static_cast<S>(0))
      : scanner(scanner),
        matrix(matrix),
        pixel_size(pixel_size),
        start_pixel(start_pixel) {}

  /// Executes Monte-Carlo system matrix generation for given detector ring
  template <class RNG,
            class AcceptanceModel,
            class ProgressCallback>
  void operator()(
      F z,                               ///< z position for calculations
      RNG& rng,                          ///< random number generator
      AcceptanceModel model,             ///< acceptance model
      int n_emissions,                   ///< number of emissions generated
      ProgressCallback& progress,        ///< progress callback
      bool o_collect_mc_matrix = true,   ///< enable matrix generation
      bool o_collect_pixel_stats = true  ///< enable pixel stats
      ) {
    if (n_emissions <= 0)
      return;

    const auto pixel_fov_radius = scanner.barrel.fov_radius() / pixel_size;
    const int pixel_fov_radius2 = pixel_fov_radius * pixel_fov_radius;

    util::random::uniform_real_distribution<F> one_dis(0, 1);
    util::random::uniform_real_distribution<F> phi_dis(0, F(M_PI));
    Distribution::SphericalDistribution<F> direction;

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
    RNG* mp_rngs = new (alloca(sizeof(RNG) * omp_get_max_threads()))
        RNG[omp_get_max_threads()];
    for (auto t = 0; t < omp_get_max_threads(); ++t) {
      mp_rngs[t].seed(rng());
    }

#pragma omp parallel for schedule(dynamic)
#endif
    // iterating only triangular matrix,
    // being upper right part or whole system matrix
    // NOTE: we must iterate pixel indices instead of x, y since we need proper
    // thread distribution when issuing on MIC
    for (int i_pixel = 0; i_pixel < matrix.total_n_pixels_in_triangle;
         ++i_pixel) {
      progress(i_pixel);
      TRY;
      auto pixel = matrix.pixel_at_index(i_pixel);

      if (pixel.x < start_pixel || pixel.y < start_pixel ||
          pixel.distance_from_origin2() > pixel_fov_radius2)
        continue;

      int pixel_hit_count = 0;
      for (auto n = 0; n < n_emissions; ++n) {

#if _OPENMP
        auto& l_rng = mp_rngs[omp_get_thread_num()];
#else
        auto& l_rng = rng;
#endif
        auto rx = (pixel.x + one_dis(l_rng)) * pixel_size;
        auto ry = (pixel.y + one_dis(l_rng)) * pixel_size;

        // ensure we are within a triangle, so we got only half hits on diagonal
        if (rx > ry)
          continue;

        auto rz = z + one_dis(l_rng) * pixel_size;
        typename ScannerClass::Response response;

        Event event(PET3D::Point<float>(rx, ry, rz), direction(l_rng));
        auto hits = scanner.detect(l_rng, model, event, response);

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

      }  // loop over emmisions from pixel

      matrix.compact_pixel_index(i_pixel);
      CATCH;
    }
  }

 private:
  const Scanner& scanner;
  Matrix& matrix;
  F pixel_size;
  F tof_step;
  S start_pixel;
};

}  // Hybrid
}  // PET3D
