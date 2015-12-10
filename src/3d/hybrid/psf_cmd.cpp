/// \page cmd_3d_hybrid_psf 3d_hybrid_psf
/// \brief 3D Longitudinal PET Point-Spread-Function calculation tool
///
/// Creates Point-Spread-Function statistics for given reconstruction / image
/// files.
///
/// Usage
/// -----
/// \verboutput 3d_hybrid_psf
///
/// \sa \ref cmd_3d_hybrid_phantom, \ref cmd_3d_hybrid_matrix

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

void print_psf(const std::string& fn,
               const VoxelMap& img,
               const F s_pixel,
               std::ostream& out);

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  PET3D::Hybrid::add_psf_options(cl);
  cl.parse_check(argc, argv);
  PET3D::Hybrid::calculate_psf_options(cl, argc);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  auto n_pixels = cl.get<int>("n-pixels");
  auto s_pixel = cl.get<double>("s-pixel");
  PixelGrid pixel_grid(n_pixels, n_pixels, s_pixel);

  int n_planes = cl.get<int>("n-planes");
  VoxelGrid grid(pixel_grid, -s_pixel * n_planes / 2, n_planes);

  std::cerr << "   voxel grid = "  // grid size:
            << grid.pixel_grid.n_columns << " x " << grid.pixel_grid.n_rows
            << " x " << grid.n_planes << std::endl;
  std::cerr << "   voxel size = " << s_pixel << std::endl;

  for (const auto& fn : cl.rest()) {
#if DEBUG
    std::cerr << "     response = " << fn << std::endl;
#endif
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
    if (!bin.is_open()) {
      throw("cannot open: " + fn);
    }
    bin >> img;
#endif
#if DEBUG
    std::cerr << "         data = " << img.data << std::endl;
#endif
    print_psf(fn, img, s_pixel, std::cout);
#if USE_MMAP
    munmap(data, data_size);
    close(fd);
#endif
  }

  CMDLINE_CATCH
}

void print_psf(const std::string& fn,
               const VoxelMap& img,
               const F s_pixel,
               std::ostream& out) {
  (void)img;
  Voxel max_voxel;
  F max;
  PET3D::Hybrid::PSF::find_max(img, max_voxel, max);
  Voxel left_above_half, right_above_half;
  PET3D::Hybrid::PSF::find_left_right_above_half(
      img, max_voxel, max, left_above_half, right_above_half);
  Vector left, right, psf;
  PET3D::Hybrid::PSF::calculate(
      img, max_voxel, max, left_above_half, right_above_half, left, right, psf);
  out << std::setw(35) << fn << '\t'         //
      << std::setw(3) << std::fixed          //
      << max_voxel.x << ' '                  //
      << max_voxel.y << ' '                  //
      << max_voxel.z << ' '                  //
      << std::setw(15) << std::fixed << max  //
      << std::setw(3) << '\t'                //
#if PRINT_ABOVE
      << left_above_half.x << ' '    //
      << left_above_half.y << ' '    //
      << left_above_half.z << '\t'   //
      << right_above_half.x << ' '   //
      << right_above_half.y << ' '   //
      << right_above_half.z << '\t'  //
#endif
      << std::setfill(' ') << std::setw(6)  //
      << std::setprecision(2)               //
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
