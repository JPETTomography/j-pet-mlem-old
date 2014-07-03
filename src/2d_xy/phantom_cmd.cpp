// PET Phantom
// Authors:
//   Piotr Bialas    <piotr.bialas@uj.edu.pl>
//   Jakub Kowal     <jakub.kowal@uj.edu.pl>
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//
// Generates phantom measurements using Monte Carlo.

#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"

#include "geometry/point.h"
#include "phantom.h"
#include "detector_ring.h"
#include "circle_detector.h"
#include "triangle_detector.h"
#include "polygonal_detector.h"
#include "model.h"
#include "util/png_writer.h"
#include "util/progress.h"

#if _OPENMP
#include <omp.h>
#endif

// all available detector shapes
typedef DetectorRing<double, int, SquareDetector<double>> SquareDetectorRing;
typedef DetectorRing<double, int, CircleDetector<double>> CircleDetectorRing;
typedef DetectorRing<double, int, TriangleDetector<double>>
    TriangleDetectorRing;
typedef DetectorRing<double, int, PolygonalDetector<6, double>>
    HexagonalDetectorRing;

template <typename DetectorRing, typename Model>
void run(cmdline::parser& cl, Model& model);

int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    cl.footer("phantom_description");

    cl.add<cmdline::path>("config",
                          'c',
                          "load config file",
                          cmdline::dontsave,
                          cmdline::path(),
                          cmdline::default_reader<cmdline::path>(),
                          cmdline::load);
#if _OPENMP
    cl.add<int>(
        "n-threads", 't', "number of OpenMP threads", cmdline::dontsave);
#endif
    cl.add<int>(
        "n-pixels", 'n', "number of pixels in one dimension", false, 256);
    cl.add<int>("m-pixel", 0, "starting pixel for partial matrix", false, 0);
    cl.add<int>("n-detectors", 'd', "number of ring detectors", false, 64);
    cl.add<int>("n-emissions",
                'e',
                "emissions",
                cmdline::optional,
                0,
                cmdline::default_reader<int>(),
                cmdline::not_from_file);
    cl.add<double>("radius", 'r', "inner detector ring radius", false);
    cl.add<double>("s-pixel", 'p', "pixel size", false);
    cl.add<double>(
        "tof-step", 'T', "TOF quantisation step for distance delta", false);
    cl.add<std::string>(
        "shape",
        'S',
        "detector (scintillator) shape (square, circle, triangle, hexagon)",
        false,
        "square",
        cmdline::oneof<std::string>("square", "circle", "triangle", "hexagon"));
    cl.add<double>("w-detector", 'w', "detector width", false);
    cl.add<double>("h-detector", 'h', "detector height", false);
    cl.add<double>("d-detector",
                   0,
                   "inscribe detector shape into circle of given diameter",
                   false);
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
                   'l',
                   "scintillator emission base length P(l)=1-e^(-1)",
                   false,
                   0.1);
    cl.add<std::mt19937::result_type>(
        "seed", 's', "random number generator seed", cmdline::dontsave);
    cl.add<cmdline::path>("output",
                          'o',
                          "output lor hits for supplied phantom",
                          cmdline::dontsave);
    cl.add("detected", 0, "collects detected emissions");

    // printing & stats params
    cl.add("verbose", 'v', "prints the iterations information on std::out");

    cl.try_parse(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    // convert obsolete acceptance option to length scale
    auto& length_scale = cl.get<double>("base-length");
    if (cl.exist("acceptance") && !cl.exist("base-length")) {
      length_scale = 1.0 / cl.get<double>("acceptance");
    }
    // FIXME: fixup for spelling mistake, present in previous versions
    auto& model = cl.get<std::string>("model");
    if (model == "scintilator") {
      model = "scintillator";
    }

    auto& n_pixels = cl.get<int>("n-pixels");
    auto& n_detectors = cl.get<int>("n-detectors");
    auto& radius = cl.get<double>("radius");
    auto& s_pixel = cl.get<double>("s-pixel");
    auto& w_detector = cl.get<double>("w-detector");
    auto& d_detector = cl.get<double>("d-detector");
    auto& shape = cl.get<std::string>("shape");
    auto verbose = cl.exist("verbose");

    if (verbose) {
      std::cerr << "assumed:" << std::endl;
    }

    // automatic radius size
    if (!cl.exist("radius")) {
      if (!cl.exist("s-pixel")) {
        radius = M_SQRT2;  // exact result
      } else {
        radius = s_pixel * n_pixels / M_SQRT2;
      }
      std::cerr << "--radius=" << radius << std::endl;
    }

    // automatic pixel size
    if (!cl.exist("s-pixel")) {
      if (!cl.exist("radius")) {
        s_pixel = 2. / n_pixels;  // exact result
      } else {
        s_pixel = M_SQRT2 * radius / n_pixels;
      }
      std::cerr << "--s-pixel=" << s_pixel << std::endl;
    }

    if (!cl.exist("w-detector")) {
      if (cl.exist("n-detectors")) {
        w_detector = 2 * M_PI * .9 * radius / n_detectors;
      } else if (cl.exist("d-detector")) {
        if (shape == "circle") {
          w_detector = d_detector;
        } else {
          auto mult = 1.;
          auto sides = 0.;
          if (shape == "triangle") {
            sides = 3.;
          } else if (shape == "square") {
            sides = 4.;
          } else if (shape == "hexagon") {
            sides = 6.;
            mult = 2.;
          } else {
            throw("cannot determine detector width for given shape");
          }
          w_detector = d_detector * std::sin(M_PI / sides) * mult;
        }
      }
      std::cerr << "--w-detector=" << w_detector << std::endl;
    }

    // automatic detector size
    // NOTE: detector height will be determined per shape
    if (!cl.exist("n-detectors")) {
      if (cl.exist("d-detector")) {
        n_detectors =
            ((int)std::floor(
                 M_PI / std::atan2(d_detector, 2 * radius + d_detector / 2)) /
             4) *
            4;
      } else {
        n_detectors =
            ((int)std::floor(M_PI / std::atan2(w_detector, 2 * radius)) / 4) *
            4;
      }
      if (!n_detectors) {
        throw("detector width is too big for given detector ring radius");
      }
      std::cerr << "--n-detectors=" << n_detectors << std::endl;
    }

    // run simmulation on given detector model & shape
    if (model == "always") {
      AlwaysAccept<> model;
      if (shape == "square") {
        run<SquareDetectorRing>(cl, model);
      } else if (shape == "circle") {
        run<CircleDetectorRing>(cl, model);
      } else if (shape == "triangle") {
        run<TriangleDetectorRing>(cl, model);
      } else if (shape == "hexagon") {
        run<HexagonalDetectorRing>(cl, model);
      }
    } else if (model == "scintillator") {
      ScintilatorAccept<> model(length_scale);
      if (shape == "square") {
        run<SquareDetectorRing>(cl, model);
      } else if (shape == "circle") {
        run<CircleDetectorRing>(cl, model);
      } else if (shape == "triangle") {
        run<TriangleDetectorRing>(cl, model);
      } else if (shape == "hexagon") {
        run<HexagonalDetectorRing>(cl, model);
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

template <typename DetectorRing, typename Model>
void run(cmdline::parser& cl, Model& model) {

  auto& n_pixels = cl.get<int>("n-pixels");
  auto& m_pixel = cl.get<int>("m-pixel");
  auto& n_detectors = cl.get<int>("n-detectors");
  auto& n_emissions = cl.get<int>("n-emissions");
  auto& radius = cl.get<double>("radius");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& w_detector = cl.get<double>("w-detector");
  auto& h_detector = cl.get<double>("h-detector");
  auto& tof_step = cl.get<double>("tof-step");
  auto verbose = cl.exist("verbose");

  PointSources<> point_sources;

  // automatic radius
  if (!cl.exist("s-pixel")) {
    if (!cl.exist("radius")) {
      s_pixel = 2. / n_pixels;  // exact result
    } else {
      s_pixel = M_SQRT2 * radius / n_pixels;
    }
    std::cerr << "--s-pixel=" << s_pixel << std::endl;
  }

  // automatic detector size
  if (!cl.exist("w-detector")) {
    w_detector = 2 * M_PI * .9 * radius / n_detectors;
    std::cerr << "--w-detector=" << w_detector << std::endl;
  }
  // NOTE: detector height will be determined per shape

  std::random_device rd;
  std::mt19937 gen(rd());
  if (cl.exist("seed")) {
    gen.seed(cl.get<std::mt19937::result_type>("seed"));
  }

  DetectorRing dr(n_detectors, radius, w_detector, h_detector);

  int n_tof_positions = 1;
  double max_bias = 0;
  if (cl.exist("tof-step") && tof_step > 0) {
    max_bias = Model::max_bias();
    n_tof_positions = dr.n_positions(tof_step, max_bias);
  }

  int* tubes = new int[n_detectors * n_detectors * n_tof_positions]();

  int pixels[n_pixels][n_pixels];
  int pixels_detected[n_pixels][n_pixels];

  for (auto i = 0; i < n_pixels; ++i)
    for (auto j = 0; j < n_pixels; ++j) {
      pixels[i][j] = 0;
      pixels_detected[i][j] = 0;
    }

  int n_emitted = 0;
  bool only_detected = false;
  if (cl.exist("detected"))
    only_detected = true;

  Phantom<> phantom;

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
        point_sources.push_back(PointSource<>(is));
      } else if (type == "ellipse") {
        phantom.push_back(EllipticalRegion<>(is));
      } else {
        std::ostringstream msg;
        msg << fn << ":" << n_line << " unhandled type of shape: " << type;
        throw(msg.str());
      }
    } while (!in.eof());
  }

  uniform_real_distribution<> one_dis(0, 1);
  uniform_real_distribution<> point_dis(-n_pixels * s_pixel / 2,
                                        +n_pixels * s_pixel / 2);
  uniform_real_distribution<> phi_dis(0, M_PI);

  Progress progress(verbose, n_emissions, only_detected ? 10000 : 1000000);

  auto fov_radius2 = dr.fov_radius() * dr.fov_radius();

  if (phantom.n_regions() > 0) {
    while (n_emitted < n_emissions) {
      progress(n_emitted);

      Point<> p(point_dis(gen), point_dis(gen));

      if (p.length2() >= fov_radius2)
        continue;

      if (phantom.test_emit(p, one_dis(gen))) {
        auto pixel = p.pixel(s_pixel, n_pixels / 2);
        // ensure we are inside pixel matrix
        if (pixel.x >= n_pixels || pixel.y >= n_pixels || pixel.x <= m_pixel ||
            pixel.y <= m_pixel)
          continue;

        typename DetectorRing::LOR lor;
        pixels[pixel.y][pixel.x]++;
        auto angle = phi_dis(gen);
        double position;
        auto hits = dr.emit_event(gen, model, p.x, p.y, angle, lor, position);
        if (hits == 2) {
          if (lor.first > lor.second)
            std::swap(lor.first, lor.second);
          int quantized_position = 0;
          if (n_tof_positions > 1) {
            quantized_position =
                dr.quantize_position(position, tof_step, n_tof_positions);
          }
          tubes[(lor.first * n_detectors + lor.second) * n_tof_positions +
                quantized_position]++;
          pixels_detected[pixel.y][pixel.x]++;
          if (only_detected)
            n_emitted++;
        }
        if (!only_detected)
          n_emitted++;
      }
    }
  }
  progress(n_emitted);

  if (point_sources.n_sources() > 0) {
    point_sources.normalize();
    n_emitted = 0;
    while (n_emitted < n_emissions) {
      progress(n_emitted);

      auto rng = one_dis(gen);
      auto p = point_sources.draw(rng);

      if (p.length2() >= fov_radius2)
        continue;

      auto pixel = p.pixel(s_pixel, n_pixels / 2);
      // ensure we are inside pixel matrix
      if (pixel.x >= n_pixels || pixel.y >= n_pixels || pixel.x <= m_pixel ||
          pixel.y <= m_pixel)
        continue;

      pixels[pixel.y][pixel.x]++;
      auto angle = phi_dis(gen);
      typename DetectorRing::LOR lor;
      double position;
      auto hits = dr.emit_event(gen, model, p.x, p.y, angle, lor, position);
      if (hits == 2) {
        if (lor.first > lor.second)
          std::swap(lor.first, lor.second);
        int quantized_position = 0;
        if (n_tof_positions > 1) {
          quantized_position =
              dr.quantize_position(position, tof_step, n_tof_positions);
        }
        tubes[(lor.first * n_detectors + lor.second) * n_tof_positions +
              quantized_position]++;
        pixels_detected[pixel.y][pixel.x]++;
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

  std::ofstream os(fn_wo_ext + ".cfg", std::ios::trunc);
  os << cl;

  png_writer pix(fn_wo_ext + ".png");
  png_writer pix_detected(fn_wo_ext + "_detected.png");

  std::ofstream pixels_text_out(fn_wo_ext + "_pixels.txt");
  std::ofstream pixels_detected_text_out(fn_wo_ext + "_detected_pixels.txt");

  int pix_max = 0;
  int pix_detected_max = 0;

  pix.write_header<>(n_pixels, n_pixels);
  pix_detected.write_header<>(n_pixels, n_pixels);

  for (auto i = 0; i < n_pixels; ++i) {
    for (auto j = 0; j < n_pixels; ++j) {
      pix_max = std::max(pix_max, pixels[i][j]);
      pix_detected_max = std::max(pix_detected_max, pixels_detected[i][j]);
    }
  }

  auto pix_gain =
      static_cast<double>(std::numeric_limits<uint8_t>::max()) / pix_max;
  for (int y = n_pixels - 1; y >= 0; --y) {
    uint8_t row[n_pixels];
    for (auto x = 0; x < n_pixels; ++x) {
      row[x] = std::numeric_limits<uint8_t>::max() - pix_gain * pixels[y][x];
    }
    pix.write_row(row);
  }

  auto pix_detected_gain =
      static_cast<double>(std::numeric_limits<uint8_t>::max()) /
      pix_detected_max;
  for (int y = n_pixels - 1; y >= 0; --y) {
    uint8_t row[n_pixels];
    for (auto x = 0; x < n_pixels; ++x) {
      row[x] = std::numeric_limits<uint8_t>::max() -
               pix_detected_gain * pixels_detected[y][x];
    }
    pix_detected.write_row(row);
  }
  for (auto i = 0; i < n_pixels; ++i) {
    for (auto j = 0; j < n_pixels; ++j) {
      pixels_text_out << pixels[i][j] << " ";
      pixels_detected_text_out << pixels_detected[i][j] << " ";
    }
    pixels_text_out << "\n";
    pixels_detected_text_out << "\n";
  }
}
