#ifdef __SSE3__
#include <xmmintrin.h>
#include <pmmintrin.h>
#endif

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"
#include "util/png_writer.h"
#include "util/backtrace.h"
#include "util/progress.h"
#include "options.h"

#include "2d/barrel/options.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/circle_detector.h"
#include "2d/barrel/square_detector.h"

#include "common/types.h"

template <class DetectorClass>
using Scanner = PET2D::Barrel::GenericScanner<DetectorClass, S>;
template <class DetectorClass>
using ScannerBuilder = PET2D::Barrel::ScannerBuilder<DetectorClass>;

using SquareScanner = Scanner<PET2D::Barrel::SquareDetector<F>>;

using Point = PET2D::Point<F>;
using Vector = PET2D::Vector<F>;

void convert_warsaw(cmdline::parser& cl);

const double cm = 0.01;
const double speed_of_light_m_per_ps = 299792458.0e-12;

int main(int argc, char* argv[]) {
  CMDLINE_TRY

#ifdef __SSE3__
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif

  cmdline::parser cl;
  cl.footer("means");
  PET2D::Barrel::add_matrix_options(cl);
  cl.add<double>("s-z", 0, "sigma z", cmdline::alwayssave, 0.01);
  cl.add<int>("only",
              '\0',
              "specifies which class (scattering) of events to include",
              false,
              1);
  cl.add<double>("length", 'l', "scintillator length", false, 1.0);
  cl.parse_check(argc, argv);
  cmdline::load_accompanying_config(cl, false);

  PET2D::Barrel::calculate_scanner_options(cl, argc);

  convert_warsaw(cl);

  CMDLINE_CATCH
}

void convert_warsaw(cmdline::parser& cl) {
  const F sigma_z = cl.get<double>("s-z");
  const F sigma_dl = cl.get<double>("s-dl");
  const F length = cl.get<double>("length");

  std::normal_distribution<double> z_error(0, sigma_z);
  std::normal_distribution<double> dl_error(0, sigma_dl);

  std::mt19937_64 rng(cl.get<long>("seed"));

  auto scanner = ScannerBuilder<SquareScanner>::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, F));

  int only = 0;
  if (cl.exist("only")) {
    only = cl.get<int>("only");
  }

  for (const auto& fn : cl.rest()) {
    std::ifstream in_means(fn);
    if (!in_means.is_open())
      throw("cannot open input file: " + cl.get<cmdline::path>("mean"));
    for (;;) {
      const int LINE_LENGTH = 256;
      char line[LINE_LENGTH];
      in_means.getline(line, LINE_LENGTH);
      if (in_means.eof())
        break;

      std::stringstream in(line);
      double x1, y1, z1, t1;
      int scatter;
      in >> x1 >> y1 >> z1 >> t1;
      double x2, y2, z2, t2;
      in >> x2 >> y2 >> z2 >> t2;
      in >> scatter;

      if (std::abs(z1 * cm) > length / 2 || std::abs(z2 * cm) > length / 2)
        continue;

      if (only == 0 || only == scatter) {

        int d1 = -1, d2 = -1;
        for (size_t i = 0; i < scanner.size(); ++i) {
          if (scanner[i].contains(Point(x1 * cm, y1 * cm), 0.0001))
            d1 = i;
          if (scanner[i].contains(Point(x2 * cm, y2 * cm), 0.0001))
            d2 = i;
        }

        if (d1 >= 0 && d2 >= 0) {

          if (d1 < d2) {
            std::swap(d1, d2);
            std::swap(z1, z2);
            std::swap(t1, t2);
          }
          double dl = (t1 - t2) * speed_of_light_m_per_ps;
          std::cout << d1 << " " << d2 << " " << z1 * cm + z_error(rng) << " "
                    << z2 * cm + z_error(rng) << " " << dl + dl_error(rng)
                    << "\n";
        } else {
          std::cerr << "Point " << Point(x1 * cm, y1 * cm) << " or "
                    << Point(x2 * cm, y2 * cm)
                    << "  outside detector - skiping\n";
        }
      }
    }
  }
}
