/// \page cmd_2d_strip_phantom 2d_strip_phantom
/// \brief 2D Strip PET phantom tool
///
/// Simulates scanner response for given virtual phantom and produces mean file
/// for \ref cmd_2d_strip_reconstruction.
///
/// Example
/// -------
///
/// 1. Make a \c playground directory and step into it
///
///        mkdir playground
///        cd playground
///
/// 2. Generate \c p_shepp_2d_strip.txt list-mode response file
///    for \c s_shepp phantom description file scaled to 30% (\c --scale)
///    and "big" barrel configuration (\c -c)
///    using additive method, 4mm pixel size (\c -p),
///    100 thousand detected emissions (\c -e together with \c --detected)
///    and be verbose (\c -v)
///
///        ../2d_strip_phantom ../phantoms/s_shepp \
///          --scale 0.3 --additive \
///          -p 0.004 \
///          -c ../config/big.cfg \
///          -e 100k --detected \
///          -o p_shepp_2d_strip.txt \
///          -v
///
/// \note Accompanying \c png files will be generated for emitted, detected and
/// naive reconstruction.
/// \note Default 1.5 cm \f$ \sigma_z \f$ and 3 cm \f$ \sigma_l \f$
/// (6 cm \f$ \sigma_{\Delta l} \f$) will be used. Different values can be
/// specified with \c --s-z and \c --s-dl respectively.
///
/// Phantom description format
/// --------------------------
///
/// 2D phantom description is a textual file, where each shape (primitive or
/// source) is specified in a separate line. There are 3 types of primivites:
/// \c ellipse, \c rectangle and \c point.
///
/// - \c ellipse (all fields after "ellipse" are float values)
///
///       ellipse x_pos y_pos a_radius b_radius rotation_deg intensity
///
/// - \c rectangle line format (all fields after "rectangle" are float values)
///
///       ellipse x_pos y_pos width height intensity
///
/// - \c point line format (all fields after "point" are float values)
///
///       point x_pos y_pos intensity
///
/// All primitives' intensities are exclusive, where the top primitive lies in
/// the 1st line, unless \c --additive option is given, which makes the
/// intensities to add (accumulate) when they overlap.
///
/// Sample phantom descriptions
/// ---------------------------
/// - Shepp like phantom
///
///   \verbinclude phantoms/s_shepp
///
/// - Small Shepp like phantom
///
///   \verbinclude phantoms/s_shepp_scaled
///
/// Output format
/// -------------
/// Phantom reponse simulation outputs triples representing \f$ \tilde{z_u},
/// \tilde{z_d}, \Delta l \f$ of a single response. Depending on given extension
/// of the output file, the output can be textual (\c txt extension) or binary
/// (any other extension). Produced output file is compatible with \ref
/// cmd_2d_strip_reconstruction.
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
#include <mutex>

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

  if (!cl.rest().size() && !cl.exist("point")) {
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

  auto n_emissions = cl.get<size_t>("n-emissions");
  auto verbose = cl.count("verbose");

  Scanner scanner(PET2D_STRIP_SCANNER_CL(cl));

  Phantom phantom(cl.get<double>("scale"), cl.exist("additive"));
  for (auto& fn : cl.rest()) {
    std::ifstream in_phantom(fn);
    phantom << in_phantom;
  }

  // direct point phantom specification
  if (cl.exist("point")) {
    auto& point_phantom_cl = cl.get<std::vector<double>>("point");
    if (point_phantom_cl.size() < 2 || point_phantom_cl.size() > 3) {
      throw("--point must specifcy 2 or 3 values x,y[,intensity]");
    }
    phantom.push_back_region(new Phantom::PointRegion(
        PET2D::Point<F>(point_phantom_cl[0], point_phantom_cl[1]),
        point_phantom_cl.size() == 3 ? point_phantom_cl[2] : 1));
  }

  phantom.calculate_cdf();

  if (verbose) {
    std::cerr << "    image size: "  //
              << scanner.n_z_pixels << " x " << scanner.n_y_pixels << std::endl
              << "    pixel size: "  //
              << scanner.pixel_width << " x " << scanner.pixel_height
              << std::endl
              << "    space size: "  //
              << scanner.size_z << " x " << scanner.size_y << std::endl
              << "space top left: "  //
              << scanner.tl_z_half_w << " x " << scanner.tl_y_half_h
              << std::endl;
  }

  MonteCarlo monte_carlo(phantom, scanner);

  RNG rng;
  Common::AlwaysAccept<F> model;

  auto output = cl.get<cmdline::path>("output");
  auto output_base_name = output.wo_ext();
  auto ext = output.ext();
  auto output_txt = ext == ".txt";
  auto no_responses = cl.exist("no-responses");
  auto full = cl.exist("full");

  std::ofstream out_wo_error, out_w_error, out_exact_events, out_full_response;
  util::obstream bin_wo_error, bin_w_error, bin_exact_events, bin_full_response;
  if (output_base_name.length() && !no_responses) {
    if (output_txt) {
      out_w_error.open(output);
      if (full) {
        out_wo_error.open(output_base_name + "_wo_error" + ext);
        out_exact_events.open(output_base_name + "_events" + ext);
        out_full_response.open(output_base_name + "_full_response" + ext);
      }
    } else {
      bin_w_error.open(output);
      if (full) {
        bin_wo_error.open(output_base_name + "_wo_error" + ext);
        bin_exact_events.open(output_base_name + "_events" + ext);
        bin_full_response.open(output_base_name + "_full_response" + ext);
      }
    }
  } else {
    no_responses = true;
  }

  auto only_detected = cl.exist("detected");
  auto n_z_pixels = cl.get<int>("n-z-pixels");
  auto n_y_pixels = cl.get<int>("n-y-pixels");
  auto s_pixel = cl.get<double>("s-pixel");

  PET2D::Point<F> image_origin(0, 0);
  if (cl.exist("image-origin")) {
    auto& image_origin_cl = cl.get<std::vector<double>>("image-origin");
    if (image_origin_cl.size() != 2) {
      throw("--image-origin must specify two values for x,y");
    }
    image_origin.x = image_origin_cl[0];
    image_origin.y = image_origin_cl[1];
  }

  PET2D::PixelGrid<F, S> pixel_grid(
      n_z_pixels,
      n_y_pixels,
      s_pixel,
      PET2D::Point<F>(-s_pixel * n_z_pixels / 2, -s_pixel * n_y_pixels / 2) +
          image_origin.as_vector());

  Image image_emitted(n_z_pixels, n_y_pixels);
  Image image_detected_exact(n_z_pixels, n_y_pixels);
  Image image_detected_w_error(n_z_pixels, n_y_pixels);

  // tangent 3D map, if tan-bins not given then we make just an 1 bin deep map,
  // but we won't write to it
  auto tan_bins = cl.get<int>("tan-bins");
  auto max_tan = (F)1.1 * scanner.scintillator_length / (2 * scanner.radius);
  PET2D::Point<F> kernel_point(0, 0);
  PET3D::VariableVoxelSizeVoxelGrid<F, S> tan_bins_grid(
      pixel_grid, -max_tan, std::max(tan_bins, 1), max_tan / tan_bins * 2);
  PET3D::VoxelMap<PET3D::Voxel<S>, Hit> tan_bins_map(
      n_z_pixels, n_y_pixels, std::max(tan_bins, 1));

  if (cl.exist("kernel-point")) {
    auto& kernel_point_cl = cl.get<std::vector<double>>("kernel-point");
    if (kernel_point_cl.size() != 2) {
      throw("--kernel-point must supply exactly 2 values for x, y");
    }
    kernel_point.x = kernel_point_cl[0];
    kernel_point.y = kernel_point_cl[1];
  }

#if _OPENMP
  std::mutex event_mutex;
#define CRITICAL std::lock_guard<std::mutex> event_lock(event_mutex);
#else
#define CRITICAL
#endif

  util::progress progress(verbose, n_emissions, 10000);
  monte_carlo(
      rng,
      model,
      n_emissions,
      [&](const Event& event) {  // emitted
        auto pixel = pixel_grid.pixel_at(event.origin);
        if (pixel_grid.contains(pixel)) {
          image_emitted.increment(pixel);
        }
      },
      [&](const Event& event, const FullResponse& full_response) {  // detected
        auto response_w_error = scanner.response_w_error(rng, full_response);
        if (!no_responses) {
          if (output_txt) {
            if (full) {
              std::ostringstream ss_wo_error, ss_w_error, ss_exact_events,
                  ss_full_response;
              ss_exact_events << event << "\n";
              ss_full_response << full_response << "\n";
              ss_wo_error << scanner.response_wo_error(full_response) << "\n";
              ss_w_error << scanner.response_w_error(rng, full_response)
                         << "\n";
              {
                CRITICAL
                out_exact_events << ss_exact_events.str();
                out_full_response << ss_full_response.str();
                out_wo_error << ss_wo_error.str();
                out_w_error << ss_w_error.str();
              }
            } else {
              std::ostringstream ss_w_error;
              ss_w_error << scanner.response_w_error(rng, full_response)
                         << "\n";
              {
                CRITICAL
                out_w_error << ss_w_error.str();
              }
            }
          } else {
            CRITICAL
            if (full) {
              bin_exact_events << event;
              bin_full_response << full_response;
              bin_wo_error << scanner.response_wo_error(full_response);
            }
            bin_w_error << scanner.response_w_error(rng, full_response);
          }
        }
        {
          auto pixel = pixel_grid.pixel_at(event.origin);
          if (pixel_grid.contains(pixel)) {
            image_detected_exact.increment(pixel);
          }
        }
        {
          auto event_w_error =
              scanner.from_projection_space_tan(response_w_error);
          auto pixel = pixel_grid.pixel_at(
              PET2D::Point<F>(event_w_error.z, event_w_error.y));
          if (pixel_grid.contains(pixel)) {
            image_detected_w_error.increment(pixel);
          }
          if (tan_bins > 0) {
            auto tan_voxel = tan_bins_grid.voxel_at(PET3D::Point<F>(
                event_w_error.z, event_w_error.y, event_w_error.tan));
            if (tan_bins_grid.contains(tan_voxel)) {
              tan_bins_map.increment(tan_voxel);
            }
          }
        }
      },
      progress,
      only_detected);
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
    util::obstream bin_detected(output_base_name + "_detected");
    util::nrrd_writer nrrd_detected(output_base_name + "_detected.nrrd",
                                    output_base_name + "_detected");
    bin_detected << image_detected_exact;
    nrrd_detected << image_detected_exact;
    util::obstream bin_emitted(output_base_name + "_emitted");
    util::nrrd_writer nrrd_emitted(output_base_name + "_emitted.nrrd",
                                   output_base_name + "_emitted");
    bin_emitted << image_emitted;
    nrrd_emitted << image_emitted;
    util::obstream bin_naive(output_base_name + "_naive");
    util::nrrd_writer nrrd_naive(output_base_name + "_naive.nrrd",
                                 output_base_name + "_naive");
    bin_naive << image_detected_w_error;
    nrrd_naive << image_detected_w_error;

    // PNG
    util::png_writer png_detected(output_base_name + "_detected.png");
    png_detected << image_detected_exact;
    util::png_writer png_emitted(output_base_name + "_emitted.png");
    png_emitted << image_emitted;
    util::png_writer png_naive(output_base_name + "_naive.png");
    png_naive << image_detected_w_error;

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
      auto pixel_size = tan_bins_grid.pixel_grid.pixel_size;
      auto bin_volume = pixel_size * pixel_size * tan_bins_grid.voxel_size;
      auto inv_kernel_point_sensitivity =
          only_detected ? (1 / scanner.sensitivity(kernel_point)) : (F)1;
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
                                       kernel_point);

            auto volume_elem =
                4 * scanner.radius * bin_volume * std::sqrt(1 + tan * tan);
            tan_kernel_map[voxel] =
                kernel_value * volume_elem * inv_kernel_point_sensitivity;
          }
        }
      }

      util::obstream bin_tan_kernel(output_base_name + "_tan_kernel");
      util::nrrd_writer nrrd_tan_kernel(output_base_name + "_tan_kernel.nrrd",
                                        output_base_name + "_tan_kernel");
      bin_tan_kernel << tan_kernel_map;
      nrrd_tan_kernel << tan_kernel_map;
    }
  }

  CMDLINE_CATCH
}
