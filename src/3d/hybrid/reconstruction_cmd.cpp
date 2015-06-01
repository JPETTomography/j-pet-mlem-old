#include <iostream>
#include <fstream>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "2d/barrel/options.h"
#include "2d/barrel/ring_scanner.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"

using FType = float;
using SType = int;
using RNGType = std::mt19937;
using Detector = PET2D::Barrel::SquareDetector<FType>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, 192, SType>;
using Point = PET2D::Point<FType>;


int main(int argc, char* argv[]) {

  cmdline::parser cl;
  PET2D::Barrel::add_matrix_options(cl);
  try {
    cl.try_parse(argc, argv);
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

  PET2D::Barrel::set_big_barrel_options(cl);
  auto scanner = PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, FType));

  auto output = cl.get<cmdline::path>("output");
  auto output_base_name = output.wo_ext();
  auto ext = output.ext();

  return 0;
}
