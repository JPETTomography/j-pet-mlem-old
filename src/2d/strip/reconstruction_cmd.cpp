/// \page cmd_2d_strip_reconstruction 2d_strip_reconstruction
/// \brief 2D Strip PET reconstruction tool
///
/// Reconstructs image using List-Mode with analytic kernel approximation from
/// physical scanner response or simulated response output from \ref
/// cmd_2d_strip_phantom.
///
/// Example
/// -------
///
/// 1. Make a \c playground directory and step into it
///
///        mkdir playground
///        cd playground
///
/// 2. Reconstruct \c p_shepp_2d_strip.txt list-mode response file generated
///    with \ref cmd_2d_strip_phantom into \c r_shepp_2d_strip.txt
///    using 20 iterations and be verbose (\c -v)
///
///        ../2d_strip_reconstruction p_shepp_2d_strip.txt \
///          -i 20 \
///          -o r_shepp_2d_strip.txt \
///          -v
///
/// \note \c p_shepp_2d_strip.cfg file will be automatically read if it exists.
/// \note Accompanying \c png files will be generated for each iteration.
/// \note Default 1.5 cm \f$ \sigma_z \f$ and 3 cm \f$ \sigma_l \f$
/// (6 cm \f$ \sigma_{\Delta l} \f$) will be used. Different values can be
/// specified with \c --s-z and \c --s-dl respectively.
///
/// \image html detector-frame-space-ellipse.pdf.png
///
/// Input format
/// ------------
/// 2D strip reconstruction works in list-mode and expects an input being a list
/// of responses - triples representing \f$ \tilde{z_u}, \tilde{z_d}, \Delta l
/// \f$. Depending on the extension of the input file, it can be textual (\c txt
/// extension) or binary (any other extension). In case of textual file, each
/// line contains single response of 3 float numbers space separated, eg.:
///
///     0.0740503 0.121193 0.221124
///     0.0447861 -0.0682753 -0.0344387
///     -0.0885979 0.230476 0.186496
///     -0.07959 0.208617 0.111321
///     0.0272275 -0.0228252 -0.392962
///     0.0561384 -0.104786 -0.358995
///     0.146026 -0.0567436 -0.128083
///     -0.197322 0.0686767 0.298622
///
/// All other extensions (or no extension) is treated as a binary files, a
/// triples of 32-bit SP floating point numbers.
///
/// In order to convert various other list-mode reponse formats into this format
/// \ref cmd_3d_tool_convert can be used.
///
/// Authors
/// -------
/// - Adam Strzelecki <adam.strzelecki@uj.edu.pl>
/// - Jakub Kowal     <jakub.kowal@uj.edu.pl>
///
/// Usage
/// -----
/// \verboutput 2d_strip_reconstruction
///
/// \sa \ref cmd_2d_strip_phantom, \ref cmd_3d_tool_convert

#include <iostream>
#include <ostream>
#include <vector>

#if SSE_FLUSH
#include <xmmintrin.h>
#endif

#include "cmdline.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/png_writer.h"
#include "util/progress.h"
#include "util/backtrace.h"

#include "options.h"
#include "response.h"
#include "reconstruction.h"
#include "gaussian_kernel.h"

#include "common/types.h"

#if HAVE_CUDA
#include "cuda/reconstruction.h"
#endif

using Kernel = PET2D::Strip::GaussianKernel<F>;
using Reconstruction = PET2D::Strip::Reconstruction<F, Kernel>;
using Scanner = PET2D::Strip::Scanner<F, S>;

void print_statistics(std::ostream& out,
                      const Reconstruction& reconstruction,
                      int n_iterations,
                      int n_blocks,
                      std::string prefix = std::string());

int main(int argc, char* argv[]) {
  CMDLINE_TRY

#if SSE_FLUSH
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  cmdline::parser cl;
  PET2D::Strip::add_reconstruction_options(cl);
  cl.add<cmdline::path>("rho", 0, "start rho (eg. existing iteration)", false);
  cl.add("rho-is-int", 0, "start rho uses integer representation (not float)");
  cl.parse_check(argc, argv);
  PET2D::Strip::calculate_scanner_options(cl, argc);

  if (!cl.rest().size()) {
    if (argc == 1) {
      std::cerr << cl.usage();
      exit(0);
    } else {
      throw("at least one responses input file expected, consult --help");
    }
  }

  cmdline::load_accompanying_config(cl);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  Scanner scanner(PET2D_STRIP_SCANNER_CL(cl));
  Reconstruction reconstruction(scanner);

  auto verbose = cl.count("verbose");
  if (verbose) {
    std::cout << "# image: " << scanner.n_y_pixels << "x" << scanner.n_z_pixels
              << std::endl;
  }

  if (cl.exist("rho")) {
    auto rho_name = cl.get<cmdline::path>("rho");
    auto rho_ext = rho_name.ext();
    auto rho_txt = rho_ext == ".txt";
    if (rho_txt) {
      std::ifstream txt(rho_name);
      txt >> reconstruction.rho;
    } else {
      util::ibstream bin(rho_name);
      if (cl.exist("rho-is-int")) {
        PET2D::PixelMap<PET2D::Pixel<S>, Hit> int_rho(
            reconstruction.rho.width, reconstruction.rho.height);
        bin >> int_rho;
        for (int i = 0; i < int_rho.size; ++i) {
          reconstruction.rho[i] = int_rho[i] + 1000;
        }
      } else {
        bin >> reconstruction.rho;
      }
    }
  }

  auto is_3d = cl.exist("3d-response");

  for (auto& fn : cl.rest()) {
    if (cmdline::path(fn).ext() == ".txt") {
#if USE_FAST_TEXT_PARSER
      reconstruction.fast_load_txt_events(fn.c_str(), is_3d);
#else
      if (is_3d) {
        throw("3D input not supported in this build");
      }
      std::ifstream in_responses(fn);
      ENSURE_IS_OPEN(in_responses, "phantom responses", fn);
      reconstruction << in_responses;
#endif
    } else {
      if (is_3d) {
        throw("3D input must have .txt extension");
      }
      util::ibstream in_responses(fn);
      ENSURE_IS_OPEN(in_responses, "phantom responses", fn);
      reconstruction << in_responses;
    }
    if (verbose) {
      std::cerr << "# read " << reconstruction.responses.size()
                << " responsess from " << fn << std::endl;
    }

    auto dl_is_time = cl.exist("dl-is-time");
    if (cl.exist("scale-length")) {
      auto scale = cl.get<double>("scale-length");
      for (auto& response : reconstruction.responses) {
        response.z_u *= scale;
        response.z_d *= scale;
        if (!dl_is_time)
          response.dl *= scale;
      }
    }
    if (dl_is_time) {
      auto scale = cl.get<double>("scale-time");
      auto speed_of_light = cl.get<double>("speed-of-light");
      for (auto& response : reconstruction.responses) {
        response.dl = response.dl * speed_of_light * scale;
      }
    }
  }

  auto n_blocks = cl.get<int>("blocks");
  auto n_iterations_in_block = cl.get<int>("iterations");
  auto n_iterations = n_blocks * n_iterations_in_block;
  auto output_name = cl.get<cmdline::path>("output");
  auto output_base_name = output_name.wo_ext();
  auto output_ext = output_name.ext();
  auto output_txt = output_ext == ".txt";

  util::progress progress(verbose, n_iterations, 1);

#if HAVE_CUDA
  if (cl.exist("gpu")) {
    PET2D::Strip::GPU::Reconstruction::run<F>(
        scanner,
        reconstruction.responses.data(),
        reconstruction.responses.size(),
        reconstruction.rho,
        n_blocks,
        n_iterations_in_block,
        [&](int iteration,
            const PET2D::Strip::GPU::Reconstruction::Output& output) {
          if (!output_base_name.length())
            return;
          auto fn = iteration >= 0
                        ? output_base_name.add_index(iteration, n_iterations)
                        : output_base_name + "_sensitivity";
          util::png_writer png(fn + ".png");
          png << output;
          if (output_txt) {
            std::ofstream txt(fn + ".txt");
            txt << output;
          } else if (output_ext != ".png") {
            util::obstream bin(fn + output_ext);
            bin << output;
          }
          util::nrrd_writer nrrd(fn + ".nrrd", fn + output_ext, output_txt);
          nrrd << output;
        },
        [&](int completed, bool finished) { progress(completed, finished); },
        cl.get<int>("cuda-device"),
        cl.get<int>("cuda-blocks"),
        cl.get<int>("cuda-threads"),
        [=](const char* device_name) {
          if (verbose) {
            std::cerr << "# device: " << device_name << std::endl;
          }
        });
  } else
#endif
  {
    if (output_base_name.length()) {
      util::png_writer png(output_base_name + "_sensitivity.png");
      png << reconstruction.sensitivity;
    }

    for (int block = 0; block < n_blocks; block++) {
      reconstruction(
          progress, n_iterations_in_block, block * n_iterations_in_block);

      if (!output_base_name.length())
        continue;
      auto fn = output_base_name.add_index((block + 1) * n_iterations_in_block,
                                           n_iterations);

      util::png_writer png(fn + ".png", cl.get<double>("png-max"));
      png << reconstruction.rho;
      if (output_txt) {
        std::ofstream txt(fn + ".txt");
        txt.precision(12);
        txt << std::fixed;
        reconstruction.output_tuples(txt);
      } else if (output_ext != ".png") {
        util::obstream bin(fn + output_ext);
        util::nrrd_writer nrrd(fn + ".nrrd", fn + output_ext);
        bin << reconstruction.rho;
        nrrd << reconstruction.rho;
      }
    }
#if USE_STATISTICS
    if (verbose) {
      print_statistics(std::cout, reconstruction, n_iterations, n_blocks, "# ");
    }
#endif
  }

  CMDLINE_CATCH
}

#if USE_STATISTICS
void print_statistics(std::ostream& out,
                      const Reconstruction& reconstruction,
                      int n_iterations,
                      int n_blocks,
                      std::string prefix) {
  size_t iterations = n_iterations * n_blocks;
  size_t events = reconstruction.stats.n_events_processed / iterations;
  size_t pixels = reconstruction.stats.n_pixels_processed / iterations;
  size_t kernels = reconstruction.stats.n_kernel_calls / iterations;
  out << prefix << "iterations: " << iterations << " "
      << "events: " << events << " "
      << "pixels: " << pixels << " "
      << "(" << (double)pixels / events << ") "
      << "kernel calls: " << kernels << " "
      << "(" << (double)kernels / events << ")" << std::endl;

  size_t bb_width_sum = reconstruction.stats.bb_width_sum / n_iterations;
  size_t bb_height_sum = reconstruction.stats.bb_height_sum / n_iterations;
  size_t bb_width2_sum = reconstruction.stats.bb_width2_sum / n_iterations;
  size_t bb_height2_sum = reconstruction.stats.bb_height2_sum / n_iterations;
  size_t bb_width_height_sum =
      reconstruction.stats.bb_width_height_sum / n_iterations;
  double avg_width = (double)bb_width_sum / events;
  double avg_height = (double)bb_height_sum / events;
  double avg_width2 = (double)bb_width2_sum / events;
  double avg_height2 = (double)bb_height2_sum / events;
  double avg_width_height = (double)bb_width_height_sum / events;
  avg_width2 -= avg_width * avg_width;
  avg_height2 -= avg_height * avg_height;
  avg_width_height -= avg_width * avg_height;
  out << prefix << "width: " << avg_width << "(" << std::sqrt(avg_width2) << ")"
      << " height: " << avg_height << "(" << std::sqrt(avg_height2) << ")  "
      << avg_width_height << std::endl;
}
#endif
