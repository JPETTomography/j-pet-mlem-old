/// \page cmd_3d_tool_psf 3d_tool_psf
/// \brief 3D Point-Spread-Function calculation tool
///
/// Calculates full width at half maximum (FWHM) for point spread function (PSF)
/// for given reconstruction / image files.
///
/// \image html FWHM.pdf.png
///
/// Usage
/// -----
/// \verboutput 3d_tool_psf
///
/// \sa \ref cmd_3d_hybrid_reconstruction

#if _OPENMP
#include <omp.h>
#endif

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/json.h"
#include "util/backtrace.h"

#include "common/types.h"

#include "2d/geometry/pixel_grid.h"
#include "../geometry/voxel_grid.h"
#include "../geometry/voxel_map.h"
#include "../geometry/vector.h"

#include "options.h"
#include "psf.h"

using PixelGrid = PET2D::PixelGrid<F, S>;
using VoxelGrid = PET3D::VoxelGrid<F, S>;
using Voxel = PET3D::Voxel<S>;
using VoxelMap = PET3D::VoxelMap<Voxel, F>;
using Vector = PET3D::Vector<F>;

#if !_WIN32
#define USE_MMAP 1
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#else
#include "util/bstream.h"
#endif

void print_psf(const cmdline::path& fn,
               const VoxelMap& img,
               const F s_pixel,
               std::ostream& out,
               const S padding = 0);

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  PET3D::Tool::add_psf_options(cl);
  cl.parse_check(argc, argv);
  PET3D::Tool::calculate_psf_options(cl, argc);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  auto n_pixels = cl.get<int>("n-pixels");
  auto s_pixel = cl.get<double>("s-pixel");
  auto padding = cl.get<int>("padding");
  PixelGrid pixel_grid(n_pixels, n_pixels, s_pixel);

  int n_planes = cl.get<int>("n-planes");
  VoxelGrid grid(pixel_grid, -s_pixel * n_planes / 2, n_planes);

  std::cerr << "   voxel grid = "  // grid size:
            << grid.pixel_grid.n_columns << " x " << grid.pixel_grid.n_rows
            << " x " << grid.n_planes << std::endl;
  std::cerr << "   voxel size = " << s_pixel << std::endl;

  for (const auto& fn : cl.rest()) {
    if (cmdline::path(fn).ext() == ".txt") {
      throw("text files are not supported by this tool: " + fn);
    }
#if USE_MMAP
    auto fd = open(fn.c_str(), O_RDONLY);
    if (fd == -1) {
      throw("cannot open: " + fn);
    }
    const auto data_size = grid.n_voxels * sizeof(F);
    F* data = (F*)mmap(NULL, data_size, PROT_READ, MAP_SHARED, fd, 0);
    VoxelMap img(
        grid.pixel_grid.n_columns, grid.pixel_grid.n_rows, grid.n_planes, data);
#else
    VoxelMap img(
        grid.pixel_grid.n_columns, grid.pixel_grid.n_rows, grid.n_planes);
    util::ibstream bin(fn);
    ENSURE_IS_OPEN(bin, "input image", fn);
    bin >> img;
#endif
    print_psf(fn, img, s_pixel, std::cout, padding);
#if USE_MMAP
    munmap(data, data_size);
    close(fd);
#endif
  }

  CMDLINE_CATCH
}

void print_psf(const cmdline::path& fn,
               const VoxelMap& img,
               const F s_pixel,
               std::ostream& out,
               const S padding) {
  (void)img;
  Voxel max_voxel;
  F max;
  PET3D::Tool::PSF::find_max(img, max_voxel, max, padding);
  Voxel left_above_half, right_above_half;
  PET3D::Tool::PSF::find_left_right_above_half(
      img, max_voxel, max, left_above_half, right_above_half);
  Vector left, right, psf;
  PET3D::Tool::PSF::calculate(
      img, max_voxel, max, left_above_half, right_above_half, left, right, psf);
  out << std::setw(35) << fn << ' '   //
      << std::setw(3) << std::fixed   //
      << fn.scan_index() << '\t'      //
      << max_voxel.x << ' '           //
      << max_voxel.y << ' '           //
      << max_voxel.z << ' '           //
      << std::setw(15) << std::fixed  //
      << max << '\t'                  //
      << std::setw(3)                 //
#if PRINT_ABOVE
      << left_above_half.x << ' '    //
      << left_above_half.y << ' '    //
      << left_above_half.z << '\t'   //
      << right_above_half.x << ' '   //
      << right_above_half.y << ' '   //
      << right_above_half.z << '\t'  //
#endif
      << std::setfill(' ') << std::setw(7)  //
      << std::setprecision(3)               //
#if PRINT_LEFT_RIGHT
      << left.x << ' '    //
      << left.y << ' '    //
      << left.z << '\t'   //
      << right.x << ' '   //
      << right.y << ' '   //
      << right.z << "  "  //
#endif
      << psf.x << ' '                   //
      << psf.y << ' '                   //
      << psf.z << '\t'                  //
      << psf.x * s_pixel * 1000 << ' '  //
      << psf.y * s_pixel * 1000 << ' '  //
      << psf.z * s_pixel * 1000 << std::endl;
}
