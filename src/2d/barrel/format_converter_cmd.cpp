#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"
#include "reconstruction.h"
#include "util/png_writer.h"
#include "util/backtrace.h"
#include "util/progress.h"
#include "options.h"
#include "2d/barrel/scanner_builder.h"

#include "generic_scanner.h"
#include "circle_detector.h"
#include "square_detector.h"

#include "common/types.h"

template <class DetectorClass>
using Scanner = PET2D::Barrel::GenericScanner<DetectorClass, S>;
template <class DetectorClass>
using ScannerBuilder = PET2D::Barrel::ScannerBuilder<DetectorClass>;

using SquareScanner = Scanner<PET2D::Barrel::SquareDetector<F>>;

using Point = PET2D::Point<F>;
using Vector = PET2D::Vector<F>;

const double cm = 0.01;
const double speed_of_light_m_per_s = 299792458.0;

int main(int argc, char* argv[]) {
  CMDLINE_TRY

#ifdef __SSE3__
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif

  cmdline::parser cl;
  cl.footer("means");
  PET2D::Barrel::add_matrix_options(cl);
  // PET2D::Barrel::add_config_option(cl);
  cl.parse_check(argc, argv);
  cmdline::load_accompanying_config(cl, false);

  PET2D::Barrel::calculate_scanner_options(cl, argc);

  auto scanner = ScannerBuilder<SquareScanner>::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, F));

  std::vector<PET2D::Barrel::CircleDetector<F>> circles;
  for (int i = 0; i < scanner.size(); i++) {
    auto circle = scanner[i].circumscribe_circle();
    circles.push_back(circle);
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
      in >> x1 >> y1 >> z1 >> t1;
      double x2, y2, z2, t2;
      in >> x2 >> y2 >> z2 >> t2;
      int d1, d2;
      for (int i = 0; i < scanner.size(); i++) {
        if ((Point(x1*cm, y1*cm) - circles[i].center).length2() < circles[i].radius2)
          d1 = i;
        if ((Point(x2*cm, y2*cm) - circles[i].center).length2() < circles[i].radius2)
          d2 = i;
      }
      if (d1 < d2) {
        std::swap(d1, d2);
        std::swap(z1, z2);
        std::swap(t1, t2);
      }
      std::cout << d1 << " " << d2 << " " << z1 * cm << " " << z2 * cm << " "
                << (t1 - t2) * speed_of_light_m_per_s << "\n";
    }
  }

  CMDLINE_CATCH
}
