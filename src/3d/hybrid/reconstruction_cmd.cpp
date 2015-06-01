#include <iostream>
#include <fstream>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "2d/barrel/options.h"
#include "2d/barrel/ring_scanner.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/lor_info.h"

using FType = float;
using SType = int;
using RNGType = std::mt19937;
using Detector = PET2D::Barrel::SquareDetector<FType>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, 192, SType>;
using Point = PET2D::Point<FType>;

int main(int argc, char* argv[]) {

  cmdline::parser cl;
  cl.add<std::string>("lor-info", '\0', "lor-pixel information", true);
  PET2D::Barrel::add_matrix_options(cl);
  try {
    cl.try_parse(argc, argv);

    PET2D::Barrel::set_big_barrel_options(cl);
    auto scanner =
        PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
            PET2D_BARREL_SCANNER_CL(cl, FType));

    auto output = cl.get<cmdline::path>("output");
    auto output_base_name = output.wo_ext();
    auto ext = output.ext();

    auto lor_info_file_name = cl.get<std::string>("lor-info");

    std::ifstream lor_info_istream(lor_info_file_name, std::ios::binary);

    FType pixel_size = 0.005;
    if (cl.exist("s-pixel"))
      pixel_size = cl.get<double>("s-pixel");
    FType fov_radius = cl.get<double>("fov-radius");
    std::cout << "fov " << fov_radius << " size " << pixel_size << "\n";
    SType n_columns = 2 * SType(std::ceil(fov_radius / pixel_size));
    SType n_rows = n_columns;
    std::cout << "cols " << n_columns << " rows " << n_rows << "\n";
    PET2D::PixelGrid<FType, SType> grid(
        n_columns,
        n_rows,
        pixel_size,
        Point(-pixel_size * n_columns / 2, -pixel_size * n_rows / 2));

    PET2D::Barrel::LorInfo<FType, SType> lor_info(scanner.size(), grid);
    lor_info.read(lor_info_istream);
    lor_info.print(std::cout);

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

  return 0;
}
