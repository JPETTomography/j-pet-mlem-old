/// \page cmd_3d_hybrid_phantom 3d_hybrid_phantom
/// \brief 3D Longitudinal PET phantom generation tool
///
/// Simulates detector response for given virtual phantom and produces mean file
/// for \ref cmd_3d_hybrid_reconstruction.
///
/// Usage
/// -----
/// \verboutput 3d_hybrid_phantom
///
/// \sa \ref cmd_3d_hybrid_matrix, \ref cmd_3d_hybrid_reconstruction

#include <random>
#include <mutex>

#include "cmdline.h"

#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/json.h"
#include "util/backtrace.h"
#include "util/progress.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "3d/geometry/phantom.h"
#include "3d/geometry/phantom_builder.h"
#include "3d/geometry/voxel.h"
#include "3d/geometry/voxel_map.h"
#include "3d/geometry/voxel_grid.h"

#include "scanner.h"
#include "options.h"

#include "common/model.h"
#include "common/phantom_monte_carlo.h"
#include "common/types.h"

using RNG = util::random::tausworthe;
using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, S>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Phantom = PET3D::Phantom<RNG, F>;
using Allways = Common::AlwaysAccept<F>;
using Scintillator = Common::ScintillatorAccept<F>;
using Point = PET3D::Point<F>;
using Vector = PET3D::Vector<F>;
using Voxel = PET3D::Voxel<S>;
using MonteCarlo = Common::PhantomMonteCarlo<Phantom, Scanner>;
using Event = MonteCarlo::Event;
using FullResponse = MonteCarlo::FullResponse;
using Image = PET3D::VoxelMap<Voxel, Hit>;
using Grid = PET3D::VoxelGrid<F, S>;

// FIXME: I don't know what is the purpose of this, but these are unused, so
// either should be removed or applied to the code.
#if HARDCODED_VALUES
namespace {
F strip_width = F(0.005);
F strip_height = F(0.019);
F strip_distance = F(0.410);
F inner_radius = (strip_distance - strip_height) / 2;
F strip_length = F(0.300);
}
#endif

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  PET3D::Hybrid::add_phantom_options(cl);
  cl.footer("phantom_description.json ...");
  cl.parse_check(argc, argv);
  PET3D::Hybrid::calculate_phantom_options(cl, argc);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  if (!cl.rest().size()) {
    if (argc == 1) {
      std::cerr << cl.usage();
      exit(0);
    } else {
      throw("at least one input phantom description expected, consult --help");
    }
  }

  auto verbose = cl.count("verbose");
  auto output = cl.get<cmdline::path>("output");
  auto output_base_name = output.wo_ext();
  auto ext = output.ext();
  auto output_txt = ext == ".txt";
  auto full = cl.exist("full");

  Scanner scanner(
      PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
          PET3D_LONGITUDINAL_SCANNER_CL(cl, F)),
      F(cl.get<double>("length")));
  scanner.set_sigmas(cl.get<double>("s-z"), cl.get<double>("s-dl"));

  if (output_base_name.length()) {
    std::ofstream out_json(output_base_name + ".json");
    out_json << json(scanner.barrel);
  }

  std::random_device rd;
  RNG rng(rd());
  if (cl.exist("seed")) {
    rng.seed(cl.get<util::random::tausworthe::seed_type>("seed"));
  }

  Phantom::RegionPtrList regions;

  for (const auto& fn : cl.rest()) {
    std::ifstream in(fn);
    if (!in.is_open()) {
      throw("could not open file: " + cl.get<std::string>("phantoms"));
    }
    json j;
    j << in;

    if (!j.is_object()) {
      throw("no JSON object in file:" + cl.get<std::string>("phantoms"));
    }

    const json& j_phantoms = j["phantoms"];
    if (!j_phantoms.is_array()) {
      throw("phantoms array missing in JSON file: " +
            cl.get<std::string>("phantoms"));
    }

    for (const auto& j_phantom : j_phantoms) {
      auto region = PET3D::create_phantom_region_from_json<RNG, F>(j_phantom);
#if DEBUG
      std::cerr << "Adding region\n";
#endif
      regions.push_back(region);
    }
  }

  auto n_emissions = cl.get<size_t>("n-emissions");
  auto only_detected = cl.exist("detected");
  auto additive = cl.exist("additive");

  Phantom phantom(regions, additive);

  Scintillator scintillator(F(cl.get<double>("base-length")));
  MonteCarlo monte_carlo(phantom, scanner);

  std::ofstream out_wo_error, out_w_error, out_exact_events, out_full_response;
  util::obstream bin_wo_error, bin_w_error, bin_exact_events, bin_full_response;
  if (output_base_name.length()) {
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
  }

#if _OPENMP
  std::mutex event_mutex;
#endif
  auto detected_block = [&](const Event& event,
                            const FullResponse& full_response) {
    if (output_txt) {
      if (full) {
        std::ostringstream ss_wo_error, ss_w_error, ss_exact_events,
            ss_full_response;
        ss_exact_events << event << "\n";
        ss_full_response << full_response << "\n";
        ss_wo_error << scanner.response_wo_error(full_response) << "\n";
        ss_w_error << scanner.response_w_error(rng, full_response) << "\n";
        {
#if _OPENMP
          std::lock_guard<std::mutex> event_lock(event_mutex);
#endif
          out_exact_events << ss_exact_events.str();
          out_full_response << ss_full_response.str();
          out_wo_error << ss_wo_error.str();
          out_w_error << ss_w_error.str();
        }
      } else {
        std::ostringstream ss_w_error;
        ss_w_error << scanner.response_w_error(rng, full_response) << "\n";
        {
#if _OPENMP
          std::lock_guard<std::mutex> event_lock(event_mutex);
#endif
          out_w_error << ss_w_error.str();
        }
      }
    } else {
#if _OPENMP
      std::lock_guard<std::mutex> event_lock(event_mutex);
#endif
      if (full) {
        bin_exact_events << event;
        bin_full_response << full_response;
        bin_wo_error << scanner.response_wo_error(full_response);
      }
      bin_w_error << scanner.response_w_error(rng, full_response);
    }
  };

  util::progress progress(verbose, n_emissions, 1000000);
  if (cl.exist("n-pixels")) {
    auto pixel_size = cl.get<double>("s-pixel");
    auto fov_radius = cl.get<double>("fov-radius");
    auto z_left = cl.get<double>("z-left");
    auto n_planes = cl.get<int>("n-planes");
    S n_columns, n_rows;
    if (!cl.exist("n-pixels")) {
      n_columns = 2 * S(std::ceil(fov_radius / pixel_size));
    } else {
      n_columns = cl.get<int>("n-pixels");
    }
    n_rows = n_columns;
    Grid grid(Grid::PixelGrid(n_columns, n_rows, pixel_size), z_left, n_planes);
    if (verbose) {
      std::cerr << "Image output:" << std::endl;
      std::cerr << "   voxel grid = "  // grid size:
                << grid.pixel_grid.n_columns << " x " << grid.pixel_grid.n_rows
                << " x " << grid.n_planes << " / " << grid.pixel_grid.pixel_size
                << std::endl;
    }
    Image img_emitted(
        grid.pixel_grid.n_columns, grid.pixel_grid.n_rows, grid.n_planes, 0);
    Image img_detected(
        grid.pixel_grid.n_columns, grid.pixel_grid.n_rows, grid.n_planes, 0);
    monte_carlo(
        rng,
        scintillator,
        n_emissions,
        [&](const Event& event) {
          auto voxel = grid.voxel_at(event.origin);
          if (grid.contains(voxel)) {
#if _OPENMP
            __atomic_add_fetch(
                img_emitted.data + grid.index(voxel), 1, __ATOMIC_SEQ_CST);
#else
            ++img_emitted[voxel];
#endif
          }
        },
        [&](const Event& event, const FullResponse& full_response) {
          auto voxel = grid.voxel_at(event.origin);
          if (grid.contains(voxel)) {
#if _OPENMP
            __atomic_add_fetch(
                img_detected.data + grid.index(voxel), 1, __ATOMIC_SEQ_CST);
#else
            ++img_detected[voxel];
#endif
          }
          detected_block(event, full_response);
        },
        progress,
        only_detected);

    // save images
    if (output_base_name.length()) {
      auto fn = output_base_name + "_emitted";
      util::obstream bin(fn);
      util::nrrd_writer nrrd(fn + ".nrrd", fn);
      bin << img_emitted;
      nrrd << img_emitted;
    }
    if (output_base_name.length()) {
      auto fn = output_base_name + "_detected";
      util::obstream bin(fn);
      util::nrrd_writer nrrd(fn + ".nrrd", fn);
      bin << img_detected;
      nrrd << img_detected;
    }
  } else {
    monte_carlo(rng,
                scintillator,
                n_emissions,
                [](const Event&) {},
                detected_block,
                progress,
                only_detected);
  }

  CMDLINE_CATCH
}
