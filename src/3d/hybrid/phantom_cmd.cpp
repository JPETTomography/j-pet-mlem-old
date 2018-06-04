/// \page cmd_3d_hybrid_phantom 3d_hybrid_phantom
/// \brief 3D Hybrid PET phantom generation tool
///
/// Simulates detector response for given virtual phantom and produces mean file
/// for \ref cmd_3d_hybrid_reconstruction.
///
/// Example
/// -------
///
/// 1. Make a \c playground directory and step into it
///
///        mkdir playground
///        cd playground
///
/// 2. Generate \c p_shepp_3d.txt list-mode response file
///    for \c s_shepp_like_w_shell.json 3D phantom description file
///    and "big" barrel configuration (\c -c)
///    using additive method, 4mm pixel size (\c -p),
///    \f$ 256^3 \f$ image space (\c -n),
///    100 thousand detected emissions (\c -e together with \c --detected)
///    and be verbose (\c -v)
///
///        ../3d_hybrid_phantom ../phantoms/s_shepp_like_w_shell.json \
///          -p 0.004 -n 256 \
///          -c ../config/big.cfg \
///          -e 100k --detected \
///          -o p_shepp_3d.txt \
///          -v
///
/// \note Accompanying image binary and \c nrrd header files will be generated
/// for emitted and detected 3D volumetric image densities. This requires
/// specyfying pixel size with \c -p. For viewing this data see \ref
/// cmd_3d_hybrid_reconstruction 3D volumetric image viewing section.
/// \note Default 1.5 cm \f$ \sigma_z \f$ and 3 cm \f$ \sigma_l \f$
/// (6 cm \f$ \sigma_{\Delta l} \f$) will be used. Different values can be
/// specified with \c --s-z and \c --s-dl respectively.
///
/// 3D phantom description format
/// -----------------------------
///
/// 3D phantom description is a \c json file, containing an array of
/// dictionaries for each primitive (source). There are 3 basic primivites
/// specified using \c type JSON dictionary key:
///
/// - \c cylinder with specific \c radius, \c height,
///   \c angular (dictionary) distribution description
///   with \c type equal to \c spherical,
///   and common \c intensity (float)
/// - \c ellipsoid with specific \c rx, \c ry, \c rz (float) radiuses,
///   \c angular (dictionary) distribution description
///   with \c type equal to \c spherical,
///   and common \c intensity (float)
/// - \c rectangular with specific \c rx, \c ry, \c rz (float) sizes
///   and common \c intensity (float)
/// - \c point with specific \c origin (3 float array),
///   \c angular (dictionary) distribution description,
///   with \c type equal to either \c spherical or \c single-direction,
///   and common \c intensity (float)
///
/// Both \c ellipsoid and \c rectangular have zero origin and are aligned to
/// the axes. For non-zero origin and/or roatation two special wrappers must be
/// used:
/// - \c translated with specific \c displacement (3 float array)
///   and common \c phantom containing translated primitive or other wrapper
/// - \c rotated with either specific \c R (9 float array) rotation matrix
///   or \c axis (3 float array) together with \c angle (float)
///   and common \c phantom containing translated primitive or other wrapper
///
/// All primitives' intensities are exclusive, where the top primitive lies in
/// the 1st line, unless \c --additive option is given, which makes the
/// intensities to add (accumulate) when they overlap.
///
/// Sample phantom descriptions
/// ---------------------------
/// - Shepp like 3D phantom
///
///   \verbinclude phantoms/s_shepp_like_w_shell.json
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
#include "2d/strip/response.h"
#include "3d/geometry/phantom.h"
#include "3d/geometry/phantom_builder.h"
#include "3d/geometry/voxel.h"
#include "3d/geometry/voxel_map.h"
#include "3d/geometry/voxel_grid.h"

#include "scanner.h"
#include "options.h"
#include "2d/barrel/options.h"

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

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;

  cl.add<int>("n-planes", 0, "number of voxels in z direction", false, 0);
  cl.add<double>("z-left", 0, "left extent in z direction", false, 0);
  cl.add<double>("fov-radius", '\0', "radius of the Field of View", false, 1);
  cl.add<double>("length", 'l', "length of the detector", false, 2);
  cl.add<double>(
      "s-z", 0, "TOF sigma along z axis", cmdline::alwayssave, 0.015);
  cl.add<double>("z-position", 'z', "position of the z plane", false, 0);
  cl.add<int>("tan-bins",
              0,
              "Bin reconstructed events using given number of tangent bins",
              cmdline::dontsave);
  cl.add("full", 0, "Emit additional full information", cmdline::dontsave);
  cl.add("no-responses", 0, "Do not emit responses", cmdline::dontsave);

  cl.footer("phantom_description ...");

  PET2D::Barrel::add_config_option(cl);

  cl.add<int>("n-pixels", 'n', "number of pixels in one dimension", false, 256);
  cl.add<int>("m-pixel", 0, "starting pixel for partial matrix", false, 0);

  cl.add<size_t>("n-emissions",
                 'e',
                 "number of emissions",
                 cmdline::optional,
                 0,
                 cmdline::default_reader<int>(),
                 cmdline::not_from_file);
  cl.add<double>("s-pixel", 'p', "pixel size", false);
  cl.add<double>(
      "tof-step", 't', "TOF quantisation step for distance delta", false);
  cl.add<double>("s-dl", 0, "TOF sigma delta-l", cmdline::alwayssave, 0.06);
  cl.add("lm", 0, "use list-mode instead number of hits in each lor position");
  cl.add<double>(
      "scale", 0, "Scale phantom with given constant", cmdline::alwayssave, 1);
  cl.add("additive", 0, "phantom regions are additive, not disjunctive");

  cl.add<std::string>(
      "model",
      'm',
      "acceptance model",
      false,
      "scintillator",
      cmdline::oneof<std::string>("always",
                                  "scintillator",
                                  /* obsolete */ "scintilator"));
  // NOTE: this options is obsolete (use base-length instead)
  cl.add<double>("acceptance",
                 'a',
                 "acceptance probability factor",
                 cmdline::dontsave | cmdline::hidden,
                 10.);
  cl.add<double>("base-length",
                 0,
                 "scintillator emission base length P(l)=1-e^(-1)",
                 false,
                 0.1);
  cl.add<cmdline::path>(
      "output", 'o', "output lor hits for supplied phantom", cmdline::dontsave);
  cl.add("detected", 0, "collects detected emissions");

  // printing & stats params
  cl.add("verbose", 'v', "prints the iterations information on std::out");
  cl.add<util::random::tausworthe::seed_type>(
      "seed", 'S', "random number generator seed", cmdline::dontsave);

  Common::add_openmp_options(cl);

  cl.add<cmdline::path>(
      "detector-file", '\0', "detector description file", false);
  cl.add<cmdline::path>(
      "detector-file-sym", '\0', "detector symmetries description file", false);

  cl.footer("phantom_description.json ...");
  cl.parse_check(argc, argv);

  std::stringstream assumed;
  auto calculate_pixel = cl.exist("n-pixels");

  auto& n_pixels = cl.get<int>("n-pixels");
  auto& n_planes = cl.get<int>("n-planes");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& z_left = cl.get<double>("z-left");

  if (calculate_pixel) {
    if (!cl.exist("n-planes")) {
      n_planes = n_pixels;
      assumed << "--n-planes=" << n_planes << std::endl;
    }
    if (!cl.exist("z-left")) {
      z_left = -n_planes * s_pixel / 2;
      assumed << "--z-left=" << z_left << std::endl;
    }
  }

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

  std::ifstream in_dets(cl.get<cmdline::path>("detector-file"));
  if (!in_dets) {
    std::cerr << "cannot open detector description file `"
              << cl.get<cmdline::path>("detector-file") << "'\n";
  }
  auto scanner2d =
      PET2D::Barrel::ScannerBuilder<Scanner2D>::deserialize(in_dets);
  in_dets.close();

  std::ifstream in_syms(cl.get<cmdline::path>("detector-file-sym"));
  auto symmetry = PET2D::Barrel::SymmetryDescriptor<S>::deserialize(in_syms);
  scanner2d.set_symmetry_descriptor(symmetry);
  std::cout << "n_detectors "
            << " " << Scanner2D::MaxDetectors << " " << scanner2d.size()
            << std::endl;
  in_syms.close();

  Scanner scanner(scanner2d, F(cl.get<double>("length")));
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

  // read phantoms
  for (const auto& fn : cl.rest()) {
    std::ifstream in(fn);
    ENSURE_IS_OPEN(in, "phantom description", fn);
    json j;
    j << in;

    if (!j.is_object()) {
      throw("no JSON object in file:" + fn);
    }

    const json& j_phantoms = j["phantoms"];
    if (!j_phantoms.is_array()) {
      throw("phantoms array missing in JSON file: " + fn);
    }

    for (const auto& j_phantom : j_phantoms) {
      auto region = PET3D::create_phantom_region_from_json<RNG, F>(j_phantom);
      regions.push_back(region);
    }
  }

  auto n_emissions = cl.get<size_t>("n-emissions");
  auto only_detected = cl.exist("detected");
  auto additive = cl.exist("additive");
  auto no_responses = cl.exist("no-responses");

  Phantom phantom(regions, additive);

  Scintillator scintillator(F(cl.get<double>("base-length")));
  MonteCarlo monte_carlo(phantom, scanner);

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
  }

#if _OPENMP
  std::mutex event_mutex;
#define CRITICAL std::lock_guard<std::mutex> event_lock(event_mutex);
#else
#define CRITICAL
#endif

  // this block is going to be executed per each detected event
  auto detected_block = [&](const Event& event,
                            const FullResponse& full_response) {
    if (no_responses)
      return;
    if (output_txt) {
      if (full) {
        std::ostringstream ss_wo_error, ss_w_error, ss_exact_events,
            ss_full_response;
        ss_exact_events << event << "\n";
        ss_full_response << full_response << "\n";
        ss_wo_error << scanner.response_wo_error(full_response) << "\n";
        ss_w_error << scanner.response_w_error(rng, full_response) << "\n";
        {
          CRITICAL
          out_exact_events << ss_exact_events.str();
          out_full_response << ss_full_response.str();
          out_wo_error << ss_wo_error.str();
          out_w_error << ss_w_error.str();
        }
      } else {
        std::ostringstream ss_w_error;
        ss_w_error << scanner.response_w_error(rng, full_response) << "\n";
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
  };

  util::progress progress(verbose, n_emissions, 1000000);
  // image generating mode
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

    // tangent 3D map, if tan-bins not given then we make just an 1 bin deep
    // map, but we won't write to it
    auto first_radius = scanner2d.min_max_radius().first;
    auto tan_bins = cl.get<int>("tan-bins");
    auto max_tan = (F)1.1 * cl.get<double>("length") / (2 * first_radius);
    Grid::PixelGrid tan_pixel_grid(n_planes, n_columns, pixel_size);
    PET3D::VariableVoxelSizeVoxelGrid<F, S> tan_bins_grid(
        tan_pixel_grid,
        -max_tan,
        std::max(tan_bins, 1),
        max_tan / tan_bins * 2);
    PET3D::VoxelMap<PET3D::Voxel<S>, Hit> tan_bins_map(
        n_planes, n_columns, std::max(tan_bins, 1));

    // start actual simulation
    monte_carlo(
        rng,
        scintillator,
        n_emissions,
        [&](const Event& event) {
          auto voxel = grid.voxel_at(event.origin);
          if (grid.contains(voxel)) {
            img_emitted.increment(voxel);
          }
        },
        [&](const Event& event, const FullResponse& full_response) {
          auto voxel = grid.voxel_at(event.origin);
          if (grid.contains(voxel)) {
            img_detected.increment(voxel);
          }
          detected_block(event, full_response);
          if (tan_bins > 1) {
            auto response = scanner.response_w_error(rng, full_response);
            PET2D::Strip::Response<F> strip_response(
                response.z_up, response.z_dn, response.dl);
            F tan, y, z;
            strip_response.calculate_tan_y_z(first_radius, tan, y, z);
            auto tan_voxel = tan_bins_grid.voxel_at(PET3D::Point<F>(z, y, tan));
            if (tan_bins_grid.contains(tan_voxel)) {
              tan_bins_map.increment(tan_voxel);
            }
          }
        },
        progress,
        only_detected);

    // save images
    if (output_base_name.length()) {
      {
        cmdline::path fn = output_base_name + "_emitted";
        util::obstream bin(fn);
        util::nrrd_writer nrrd(fn + ".nrrd", fn.wo_path());
        bin << img_emitted;
        nrrd << img_emitted;
      }
      {
        cmdline::path fn = output_base_name + "_detected";
        util::obstream bin(fn);
        util::nrrd_writer nrrd(fn + ".nrrd", fn.wo_path());
        bin << img_detected;
        nrrd << img_detected;
      }
      if (tan_bins > 0) {
        util::obstream bin_tan_bins(output_base_name + "_tan_bins");
        util::nrrd_writer nrrd_tan_bins(output_base_name + "_tan_bins.nrrd",
                                        output_base_name + "_tan_bins");
        bin_tan_bins << tan_bins_map;
        nrrd_tan_bins << tan_bins_map;
      }
    }
  } else {
    // no image mode
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
