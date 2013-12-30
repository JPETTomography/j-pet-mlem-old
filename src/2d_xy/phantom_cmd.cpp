// PET Phantom
// Authors:
//   Piotr Bialas    <piotr.bialas@uj.edu.pl>
//   Jakub Kowal     <jakub.kowal@uj.edu.pl>
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//
// Generates phantom measurements using Monte Carlo.

#include <random>
#include <iostream>
#include <fstream>

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

    cl.add<cmdline::string>("config",
                            'c',
                            "load config file",
                            cmdline::dontsave,
                            cmdline::string(),
                            cmdline::default_reader<cmdline::string>(),
                            cmdline::load);
#if _OPENMP
    cl.add<int>(
        "n-threads", 't', "number of OpenMP threads", cmdline::dontsave);
#endif
    cl.add<int>(
        "n-pixels", 'n', "number of pixels in one dimension", false, 256);
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
    cl.add<std::string>(
        "model",
        'm',
        "acceptance model",
        false,
        "scintillator",
        cmdline::oneof<std::string>(
            "always", "scintillator", /* obsolete */ "scintilator"));
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
    cl.add<cmdline::string>("output",
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

    // run simmulation on given detector model & shape
    auto& shape = cl.get<std::string>("shape");
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
  }
  catch (cmdline::exception& ex) {
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
  }
  catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  return 1;
}

template <typename DetectorRing, typename Model>
void run(cmdline::parser& cl, Model& model) {

  auto& n_pixels = cl.get<int>("n-pixels");
  auto& n_detectors = cl.get<int>("n-detectors");
  auto& n_emissions = cl.get<int>("n-emissions");
  auto& radius = cl.get<double>("radius");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& w_detector = cl.get<double>("w-detector");
  auto& h_detector = cl.get<double>("h-detector");
  auto& tof_step = cl.get<double>("tof-step");

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
  if (!cl.exist("h-detector")) {
    h_detector = w_detector;
    std::cerr << "--h-detector=" << h_detector << std::endl;
  }

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

  uniform_real_distribution<> one_dis(0., 1.);
  uniform_real_distribution<> fov_dis(-dr.fov_radius(), dr.fov_radius());
  uniform_real_distribution<> phi_dis(0., M_PI);

  double fov_r2 = dr.fov_radius() * dr.fov_radius();

  Phantom phantom;

  for (auto fn = cl.rest().begin(); fn != cl.rest().end(); ++fn) {
    std::ifstream in(*fn);
    if (!in.is_open()) {
      throw("cannot open input file: " + *fn);
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
        double x, y, intensity;
        is >> x >> y >> intensity;
        point_sources.add(x, y, intensity);
      } else if (type == "ellipse") {
        double x, y, a, b, angle, acceptance;
        is >> x >> y >> a >> b >> angle >> acceptance;

        phantom.add_region(x, y, a, b, angle, acceptance);
      } else {
        std::ostringstream msg;
        msg << *fn << ":" << n_line << " unhandled type of shape: " << type;
        throw(msg.str());
      }
    } while (!in.eof());
  }

  if (phantom.n_regions() > 0) {
    while (n_emitted < n_emissions) {

      double x = fov_dis(gen);
      double y = fov_dis(gen);
#if DEBUG
      std::cerr << n_emitted << " (" << x << "," << y << ")" << std::endl;
#endif
      if (x * x + y * y < fov_r2) {
        if (phantom.test_emit(x, y, one_dis(gen))) {
          auto pixel = dr.pixel(x, y, s_pixel);
          typename DetectorRing::LOR lor;
          pixels[pixel.y][pixel.x]++;
          double angle = phi_dis(gen);
          double position;
          auto hits = dr.emit_event(gen, model, x, y, angle, lor, position);
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
  }

  if (point_sources.n_sources() > 0) {
    point_sources.normalize();
    n_emitted = 0;
    while (n_emitted < n_emissions) {
      double rng = one_dis(gen);
      Point<> p = point_sources.draw(rng);

      auto pixel = dr.pixel(p.x, p.y, s_pixel);

      pixels[pixel.y][pixel.x]++;
      double angle = phi_dis(gen);
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

  auto fn = cl.get<cmdline::string>("output");
  auto fn_sep = fn.find_last_of("\\/");
  auto fn_ext = fn.find_last_of(".");
  auto fn_wo_ext =
      fn.substr(0,
                fn_ext != std::string::npos &&
                        (fn_sep == std::string::npos || fn_sep < fn_ext)
                    ? fn_ext
                    : std::string::npos);

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