/// \page cmd_2d_strip_phantom 2d_strip_phantom
/// \brief 2D Strip PET phantom tool
///
/// Simulates scanner response for given virtual phantom and produces mean file
/// for \ref cmd_2d_strip_reconstruction.
///
/// Example phantom descriptions
/// ----------------------------
/// - Shepp like phantom
///
///   \verbinclude phantoms/s_shepp
///
/// - Small Shepp like phantom
///
///   \verbinclude phantoms/s_shepp_scaled
///
/// Authors
/// -------
/// - Adam Strzelecki <adam.strzelecki@uj.edu.pl>
/// - Jakub Kowal     <jakub.kowal@uj.edu.pl>
///
/// Usage
/// -----
/// \verboutput 2d_strip_phantom
///
/// \sa \ref cmd_2d_strip_reconstruction

#include <iostream>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <random>

#if SSE_FLUSH
#include <xmmintrin.h>
#endif

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/random.h"
#include "util/backtrace.h"
#include "util/progress.h"
#include "util/png_writer.h"
#include "options.h"

#include "2d/geometry/phantom.h"
#include "2d/geometry/pixel_grid.h"
#include "2d/geometry/pixel_map.h"
#include "2d/strip/scanner.h"

#include "3d/geometry/voxel_grid.h"
#include "3d/geometry/voxel_map.h"

#include "common/model.h"
#include "common/phantom_monte_carlo.h"
#include "common/types.h"

#include "analytic_kernel.h"

#if _OPENMP
#include <omp.h>
#endif

using RNG = util::random::tausworthe;
using Scanner = PET2D::Strip::Scanner<F, S>;
using Phantom = PET2D::Phantom<RNG, F>;
using Ellipse = PET2D::Ellipse<F>;
using Image = PET2D::PixelMap<PET2D::Pixel<S>, Hit>;
using MonteCarlo = Common::PhantomMonteCarlo<Phantom, Scanner>;
using Event = MonteCarlo::Event;
using FullResponse = MonteCarlo::FullResponse;

int main(int argc, char* argv[]) {
  CMDLINE_TRY

#if SSE_FLUSH
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  cmdline::parser cl;
  PET2D::Strip::add_phantom_options(cl);
  cl.parse_check(argc, argv);
  PET2D::Strip::calculate_scanner_options(cl, argc);

  if (!cl.rest().size()) {
    if (argc == 1) {
      std::cerr << cl.usage();
      exit(0);
    } else {
      throw("at least one input phantom description expected, consult --help");
    }
  }

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  auto n_emissions = cl.get<int>("n-emissions");
  auto verbose = cl.count("verbose");

  Scanner scanner(PET2D_STRIP_SCANNER_CL(cl));

  if (verbose) {
    std::cerr << "size: " << scanner.n_z_pixels << "x" << scanner.n_y_pixels
              << std::endl;
  }

  Phantom phantom(cl.get<double>("scale"));
  for (auto& fn : cl.rest()) {
    std::ifstream in_phantom(fn);
    phantom << in_phantom;
  }

  phantom.calculate_cdf();

  if (verbose) {
    std::cerr << "scanner: " << scanner.size_y << ' ' << scanner.tl_y_half_h
              << std::endl;
  }

  MonteCarlo monte_carlo(phantom, scanner);

  RNG rng;
  Common::AlwaysAccept<F> model;

  auto output = cl.get<cmdline::path>("output");
  auto output_base_name = output.wo_ext();
  auto ext = output.ext();

  std::ofstream out_wo_error, out_w_error, out_exact_events, out_full_response;
  if (output_base_name.length()) {
    out_wo_error.open(output_base_name + "_wo_error" + ext);
    out_w_error.open(output);
    out_exact_events.open(output_base_name + "_events" + ext);
    out_full_response.open(output_base_name + "_full_response" + ext);
  }

  auto n_z_pixels = cl.get<int>("n-z-pixels");
  auto n_y_pixels = cl.get<int>("n-y-pixels");
  auto s_pixel = cl.get<double>("s-pixel");
  PET2D::PixelGrid<F, S> pixel_grid(
      n_z_pixels,
      n_y_pixels,
      s_pixel,
      PET2D::Point<F>(-s_pixel * n_z_pixels / 2, -s_pixel * n_y_pixels / 2));

  Image image_emitted(n_z_pixels, n_y_pixels);
  Image image_detected_exact(n_z_pixels, n_y_pixels);
  Image image_detected_w_error(n_z_pixels, n_y_pixels);

  // tangent 3D map, if tan-bins not given then we make just an 1 bin deep map,
  // but we won't write to it
  auto tan_bins = cl.get<int>("tan-bins");
  auto max_tan = scanner.scintillator_length / (2 * scanner.radius);
  PET3D::VoxelGrid<F, S> tan_bins_grid(
      pixel_grid, -max_tan, tan_bins > 0 ? tan_bins : 1);
  PET3D::VoxelMap<PET3D::Voxel<S>, Hit> tan_bins_map(
      n_z_pixels, n_y_pixels, tan_bins > 0 ? tan_bins : 1);

  util::progress progress(verbose, n_emissions, 10000);
  monte_carlo(
      rng,
      model,
      n_emissions,
      [&](const Event& event) {
        auto pixel = pixel_grid.pixel_at(event.origin);
        if (pixel_grid.contains(pixel)) {
          image_emitted[pixel]++;
        }
      },
      [&](const Event& event, const FullResponse& full_response) {
        out_exact_events << event << "\n";
        out_full_response << full_response << "\n";
        out_wo_error << scanner.response_wo_error(full_response) << "\n";
        auto response_w_error = scanner.response_w_error(rng, full_response);
        out_w_error << response_w_error << "\n";
        {
          auto pixel = pixel_grid.pixel_at(event.origin);
          if (pixel_grid.contains(pixel)) {
            image_detected_exact[pixel]++;
          }
        }
        {
          auto event_w_error =
              scanner.from_projection_space_tan(response_w_error);
          auto pixel = pixel_grid.pixel_at(
              PET2D::Point<F>(event_w_error.z, event_w_error.y));
          if (pixel_grid.contains(pixel)) {
            image_detected_w_error[pixel]++;
          }
          if (tan_bins > 0) {
            auto tan_voxel = tan_bins_grid.voxel_at(PET3D::Point<F>(
                event_w_error.z, event_w_error.y, event_w_error.tan));
            if (tan_bins_grid.contains(tan_voxel)) {
              tan_bins_map[tan_voxel]++;
            }
          }
        }
      },
      progress);
  if (verbose) {
    std::cerr << " emitted: " << monte_carlo.n_events_emitted() << " events"
              << std::endl
              << "detected: " << monte_carlo.n_events_detected() << " events"
              << std::endl;
  }

  if (output_base_name.length()) {
    std::ofstream cfg(output_base_name + ".cfg");
    cfg << cl;

    // RAW + NRRD
    util::obstream bin_wo_error(output_base_name + "_wo_error");
    util::nrrd_writer nrrd_wo_error(output_base_name + "_wo_error.nrrd",
                                    output_base_name + "_wo_error");
    bin_wo_error << image_detected_exact;
    nrrd_wo_error << image_detected_exact;
    util::obstream bin_emitted(output_base_name + "_emitted");
    util::nrrd_writer nrrd_emitted(output_base_name + "_emitted.nrrd",
                                   output_base_name + "_emitted");
    bin_emitted << image_emitted;
    nrrd_emitted << image_emitted;
    util::obstream bin_w_error(output_base_name + "_w_error");
    util::nrrd_writer nrrd_w_error(output_base_name + "_w_error.nrrd",
                                   output_base_name + "_w_error");
    bin_w_error << image_detected_w_error;
    nrrd_w_error << image_detected_w_error;

    if (tan_bins > 0) {
      util::obstream bin_tan_bins(output_base_name + "_tan_bins");
      util::nrrd_writer nrrd_tan_bins(output_base_name + "_tan_bins.nrrd",
                                      output_base_name + "_tan_bins");
      bin_tan_bins << tan_bins_map;
      nrrd_tan_bins << tan_bins_map;

      if (verbose) {
        std::cerr << "generating kernel map..." << std::endl;
      }
      PET3D::VoxelMap<PET3D::Voxel<S>, F> tan_kernel_map(
          n_z_pixels, n_y_pixels, tan_bins);
      for (S z = 0; z < tan_kernel_map.depth; ++z) {
        auto tan = tan_bins_grid.center_at(PET3D::Voxel<S>(0, 0, z)).z;
        auto sec = compat::sqrt(1 + tan * tan);
        for (S y = 0; y < tan_kernel_map.height; ++y) {
          for (S x = 0; x < tan_kernel_map.width; ++x) {
            PET3D::Voxel<S> voxel(x, y, z);
            auto point = tan_bins_grid.center_at(voxel);
            PET2D::Strip::AnalyticKernel<F> kernel(scanner.sigma_z,
                                                   scanner.sigma_dl);
            auto kernel_value = kernel(PET2D::Point<F>(point.x, point.y),
                                       tan,
                                       sec,
                                       scanner.radius,
                                       PET2D::Point<F>(0, 0));
            tan_kernel_map[voxel] = kernel_value;
          }
        }
      }

      util::obstream bin_tan_kernel(output_base_name + "_tan_kernel");
      util::nrrd_writer nrrd_tan_kernel(output_base_name + "_tan_kernel.nrrd",
                                        output_base_name + "_tan_kernel");
      bin_tan_kernel << tan_kernel_map;
      nrrd_tan_kernel << tan_kernel_map;
    }

    // PNG
    util::png_writer png_wo_error(output_base_name + "_wo_error.png");
    png_wo_error << image_detected_exact;
    util::png_writer png_emitted(output_base_name + "_emitted.png");
    png_emitted << image_emitted;
    util::png_writer png_w_error(output_base_name + "_w_error.png");
    png_w_error << image_detected_w_error;
  }

  CMDLINE_CATCH
}
