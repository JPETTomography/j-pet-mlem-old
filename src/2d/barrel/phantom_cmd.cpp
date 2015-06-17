/// \page cmd_2d_barrel_phantom 2d_barrel_phantom
/// \brief 2D Barrel PET phantom generation tool
///
/// Simulates detector response for given virtual phantom and produces mean file
/// for \ref cmd_2d_barrel_reconstruction.
///
/// Example phantom descriptions
/// ----------------------------
/// - Shepp like phantom
///
///   \verbinclude phantom/s_shepp
///
/// - Small Shepp like phantom
///
///   \verbinclude phantom/s_shepp_small
///
/// Authors
/// -------
/// - Piotr Bialas    <piotr.bialas@uj.edu.pl>
/// - Jakub Kowal     <jakub.kowal@uj.edu.pl>
/// - Adam Strzelecki <adam.strzelecki@uj.edu.pl>
///
/// Usage
/// -----
/// \verbinclude src/2d/barrel/phantom_cmd.txt
///
/// \sa \ref cmd_2d_barrel_matrix, \ref cmd_2d_barrel_reconstruction

#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"

#include "2d/geometry/point.h"
#include "2d/barrel/scanner_builder.h"
#include "phantom.h"
#include "ring_scanner.h"
#include "generic_scanner.h"
#include "circle_detector.h"
#include "triangle_detector.h"
#include "polygonal_detector.h"
#include "model.h"
#include "util/png_writer.h"
#include "util/progress.h"
#include "options.h"

#if _OPENMP
#include <omp.h>
#endif

using namespace PET2D;
using namespace PET2D::Barrel;

template <typename DetectorType>
using DetectorModel = GenericScanner<DetectorType, MAX_DETECTORS, short>;
// using DetectorModel = Scanner<DetectorType>;

// all available detector shapes
using SquareScanner = DetectorModel<SquareDetector<float>>;
using CircleScanner = DetectorModel<CircleDetector<float>>;
using TriangleScanner = DetectorModel<TriangleDetector<float>>;
using HexagonalScanner = DetectorModel<PolygonalDetector<6, float>>;
using PixelType = PET2D::Pixel<short>;

template <typename Scanner, typename Model>
void run(cmdline::parser& cl, Model& model);


int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    add_phantom_options(cl);
    cl.add("small", 0, "small barrel", false);
    cl.add("big", 0, "big barrel", false);
    cl.try_parse(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    if (cl.exist("big"))
      set_big_barrel_options(cl);
    calculate_scanner_options(cl);

    const auto& shape = cl.get<std::string>("shape");
    const auto& model_name = cl.get<std::string>("model");
    const auto& length_scale = cl.get<double>("base-length");

    // run simmulation on given detector model & shape
    if (model_name == "always") {
      AlwaysAccept<> model;
      if (shape == "square") {
        run<SquareScanner>(cl, model);
      } else if (shape == "circle") {
        run<CircleScanner>(cl, model);
      } else if (shape == "triangle") {
        run<TriangleScanner>(cl, model);
      } else if (shape == "hexagon") {
        run<HexagonalScanner>(cl, model);
      }
    } else if (model_name == "scintillator") {
      ScintillatorAccept<> model(length_scale);
      if (shape == "square") {
        run<SquareScanner>(cl, model);
      } else if (shape == "circle") {
        run<CircleScanner>(cl, model);
      } else if (shape == "triangle") {
        run<TriangleScanner>(cl, model);
      } else if (shape == "hexagon") {
        run<HexagonalScanner>(cl, model);
      }
    }

    return 0;
  } catch (cmdline::exception& ex) {
    if (ex.help()) {
      std::cerr << ex.usage();
    }
    for (auto& msg : ex.errors()) {
      auto name = ex.name();
      if (name) {
        std::cerr << "error at " << name << ": " << msg << std::endl;
      } else {
        std::cerr << "error: " << msg << std::endl;
      }
    }
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  return 1;
}

template <typename Detector, typename Model>
void run(cmdline::parser& cl, Model& model) {

  auto& n_pixels = cl.get<int>("n-pixels");
  auto& m_pixel = cl.get<int>("m-pixel");
  // auto& n_detectors = cl.get<int>("n-detectors");
  auto& n_emissions = cl.get<int>("n-emissions");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& tof_step = cl.get<double>("tof-step");
  auto verbose = cl.exist("verbose");

  // NOTE: detector height will be determined per shape

  std::random_device rd;
  std::mt19937 gen(rd());
  if (cl.exist("seed")) {
    gen.seed(cl.get<std::mt19937::result_type>("seed"));
  }

  Detector dr = ScannerBuilder<Detector>::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, typename Detector::F));

  auto n_detectors = dr.size();
  int n_tof_positions = 1;
  double max_bias = 0;
  if (cl.exist("tof-step") && tof_step > 0) {
    max_bias = Model::max_bias();
    n_tof_positions = dr.n_tof_positions(tof_step, max_bias);
  }

  int* tubes = new int[n_detectors * n_detectors * n_tof_positions]();

  int pixels_size = n_pixels * n_pixels * sizeof(int);
  int* pixels = (int*)alloca(pixels_size);
  int* pixels_detected = (int*)alloca(pixels_size);
  memset(pixels, 0, pixels_size);
  memset(pixels_detected, 0, pixels_size);

  int n_emitted = 0;
  bool only_detected = false;
  if (cl.exist("detected"))
    only_detected = true;

  PointPhantom<float> point_phantom;
  Phantom<float> phantom;

  for (auto& fn : cl.rest()) {
    std::ifstream in(fn);
    if (!in.is_open()) {
      throw("cannot open input file: " + fn);
    }

    int n_line = 0;
    do {
      std::string line;
      std::getline(in, line);
      ++n_line;

      if (!line.size() || line[0] == '#')
        continue;

      std::istringstream is(line);

      std::string type;
      is >> type;
      if (type == "point") {
        point_phantom.emplace_back(is);
      } else if (type == "ellipse") {
        float x, y, a, b, angle, intensity;
        // is >> x >> y >> a >> b >> angle>>intensity;
        x = util::read<float>(is);
        y = util::read<float>(is);
        a = util::read<float>(is);
        b = util::read<float>(is);
        angle = util::read<float>(is);
        intensity = util::read<float>(is);
        phantom.emplace_back(x, y, a, b, angle, intensity);
      } else {
        std::ostringstream msg;
        msg << fn << ":" << n_line << " unhandled type of shape: " << type;
        throw(msg.str());
      }
    } while (!in.eof());
  }

  util::random::uniform_real_distribution<> one_dis(0, 1);
  util::random::uniform_real_distribution<> point_dis(-n_pixels * s_pixel / 2,
                                                      +n_pixels * s_pixel / 2);
  util::random::uniform_real_distribution<> phi_dis(0, M_PI);

  util::progress progress(
      verbose, n_emissions, only_detected ? 10000 : 1000000);

  auto fov_radius2 = dr.fov_radius() * dr.fov_radius();

  if (phantom.n_regions() > 0) {
    while (n_emitted < n_emissions) {

      progress(n_emitted);

      Point<float> p(point_dis(gen), point_dis(gen));

      if (p.distance_from_origin2() >= fov_radius2)
        continue;

      if (phantom.test_emit(p, one_dis(gen))) {

        auto pixel = p.pixel<PixelType>(s_pixel, n_pixels / 2);
        if (pixel.x >= n_pixels || pixel.y >= n_pixels || pixel.x <= m_pixel ||
            pixel.y <= m_pixel)
          continue;

        // typename Detector::LOR lor;
        pixels[pixel.y * n_pixels + pixel.x]++;
        auto angle = phi_dis(gen);
        // double position;
        typename Detector::Event event(p, angle);
        typename Detector::Response response;
        auto hits = dr.detect(gen, model, event, response);
        if (hits == 2) {
          if (response.lor.first > response.lor.second)
            std::swap(response.lor.first, response.lor.second);
          int quantized_position = 0;
          if (n_tof_positions > 1) {
            quantized_position = Detector::quantize_tof_position(
                response.dl, tof_step, n_tof_positions);
          }
          tubes[(response.lor.first * n_detectors + response.lor.second) *
                    n_tof_positions +
                quantized_position]++;
          pixels_detected[pixel.y * n_pixels + pixel.x]++;
          if (only_detected)
            n_emitted++;
        }
        if (!only_detected)
          n_emitted++;
      }
    }
  }
  progress(n_emitted);

  if (point_phantom.n_sources() > 0) {
    point_phantom.normalize();
    n_emitted = 0;
    while (n_emitted < n_emissions) {
      progress(n_emitted);

      auto rng = one_dis(gen);
      auto p = point_phantom.draw(rng);

      if (p.distance_from_origin2() >= fov_radius2)
        continue;

      auto pixel = p.pixel<PixelType>(s_pixel, n_pixels / 2);
      // ensure we are inside pixel matrix
      if (pixel.x >= n_pixels || pixel.y >= n_pixels || pixel.x <= m_pixel ||
          pixel.y <= m_pixel)
        continue;

      pixels[pixel.y * n_pixels + pixel.x]++;
      auto angle = phi_dis(gen);
      // typename Detector::LOR lor;
      // double position;
      typename Detector::Event event(p, angle);
      typename Detector::Response response;
      auto hits = dr.detect(gen, model, event, response);
      if (hits == 2) {
        if (response.lor.first > response.lor.second)
          std::swap(response.lor.first, response.lor.second);
        int quantized_position = 0;
        if (n_tof_positions > 1) {
          quantized_position = Detector::quantize_tof_position(
              response.dl, tof_step, n_tof_positions);
        }
        tubes[(response.lor.first * n_detectors + response.lor.second) *
                  n_tof_positions +
              quantized_position]++;
        pixels_detected[pixel.y * n_pixels + pixel.x]++;
        if (only_detected)
          n_emitted++;
      }
      if (!only_detected)
        n_emitted++;
    }
  }

  auto fn = cl.get<cmdline::path>("output");
  auto fn_wo_ext = fn.wo_ext();
  std::ofstream n_stream(fn);

  if (n_tof_positions <= 1) {
    for (int i = 0; i < n_detectors; i++) {
      for (int j = i + 1; j < n_detectors; j++) {
        auto hits = tubes[i * n_detectors + j];
        if (hits > 0)
          n_stream << i << " " << j << "  " << hits << std::endl;
      }
    }
  } else {  // TOF
    for (int i = 0; i < n_detectors; i++) {
      for (int j = i + 1; j < n_detectors; j++) {
        for (int p = 0; p < n_tof_positions; p++) {
          auto hits = tubes[(i * n_detectors + j) * n_tof_positions + p];
          if (hits > 0)
            n_stream << i << " " << j << " " << p << "  " << hits << std::endl;
        }
      }
    }
  }

  delete[] tubes;

  std::ofstream json(fn_wo_ext + ".json", std::ios::trunc);
  dr.to_json(json);

  std::ofstream os(fn_wo_ext + ".cfg", std::ios::trunc);
  os << cl;

  util::png_writer pix(fn_wo_ext + ".png");
  util::png_writer pix_detected(fn_wo_ext + "_detected.png");

  std::ofstream pixels_text_out(fn_wo_ext + "_pixels.txt");
  std::ofstream pixels_detected_text_out(fn_wo_ext + "_detected_pixels.txt");

  int pix_max = 0;
  int pix_detected_max = 0;

  pix.write_header<>(n_pixels, n_pixels);
  pix_detected.write_header<>(n_pixels, n_pixels);

  for (auto p = 0; p < n_pixels * n_pixels; ++p) {
    pix_max = std::max(pix_max, pixels[p]);
    pix_detected_max = std::max(pix_detected_max, pixels_detected[p]);
  }

  uint8_t* row = (uint8_t*)alloca(n_pixels);

  auto pix_gain =
      static_cast<double>(std::numeric_limits<uint8_t>::max()) / pix_max;
  for (int y = n_pixels - 1; y >= 0; --y) {
    for (auto x = 0; x < n_pixels; ++x) {
      row[x] = std::numeric_limits<uint8_t>::max() -
               pix_gain * pixels[y * n_pixels + x];
    }
    pix.write_row(row);
  }

  auto pix_detected_gain =
      static_cast<double>(std::numeric_limits<uint8_t>::max()) /
      pix_detected_max;
  for (int y = n_pixels - 1; y >= 0; --y) {
    for (auto x = 0; x < n_pixels; ++x) {
      row[x] = std::numeric_limits<uint8_t>::max() -
               pix_detected_gain * pixels_detected[y * n_pixels + x];
    }
    pix_detected.write_row(row);
  }
  for (auto y = 0; y < n_pixels; ++y) {
    for (auto x = 0; x < n_pixels; ++x) {
      pixels_text_out << pixels[y * n_pixels + x] << " ";
      pixels_detected_text_out << pixels_detected[y * n_pixels + x] << " ";
    }
    pixels_text_out << "\n";
    pixels_detected_text_out << "\n";
  }
}
