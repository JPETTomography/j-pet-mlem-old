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

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
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

inline void find_max(const VoxelMap& img, Voxel& max_voxel, F& max) {
  Voxel thread_max_voxels[omp_get_max_threads()];
  F thread_maxes[omp_get_max_threads()];
  std::memset(thread_maxes, 0, sizeof(F) * omp_get_max_threads());
#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (S z = 0; z < img.depth; ++z) {
    Voxel l_max_voxel = thread_max_voxels[omp_get_thread_num()];
    auto l_max = thread_maxes[omp_get_thread_num()];
    for (S y = 0; y < img.height; ++y) {
      for (S x = 0; x < img.width; ++x) {
        const Voxel voxel(x, y, z);
        const auto value = img[voxel];
        if (value > l_max) {
          l_max_voxel = voxel;
          l_max = value;
        }
      }
    }
    thread_max_voxels[omp_get_thread_num()] = l_max_voxel;
    thread_maxes[omp_get_thread_num()] = l_max;
  }
  max = 0;
  for (int t = 0; t < omp_get_max_threads(); ++t) {
    if (thread_maxes[t] > max) {
      max = thread_maxes[t];
      max_voxel = thread_max_voxels[t];
    }
  }
}

inline void find_left_right_above_half(const VoxelMap& img,
                                       const Voxel max_voxel,
                                       const F max,
                                       Voxel& left_above_half,
                                       Voxel& right_above_half) {
  const auto half_max = max / 2;
#if _OPENMP
#pragma omp task shared(img, left_above_half, right_above_half)
#endif
  {
    for (int x = 0; x <= max_voxel.x; ++x) {
      Voxel voxel(x, max_voxel.y, max_voxel.z);
      if (img[voxel] >= half_max) {
        left_above_half.x = x;
        break;
      }
    }
    for (int x = img.width - 1; x >= max_voxel.x; --x) {
      Voxel voxel(x, max_voxel.y, max_voxel.z);
      if (img[voxel] >= half_max) {
        right_above_half.x = x;
        break;
      }
    }
  }

#if _OPENMP
#pragma omp task shared(img, left_above_half, right_above_half)
#endif
  {
    for (int y = 0; y <= max_voxel.y; ++y) {
      Voxel voxel(max_voxel.x, y, max_voxel.z);
      if (img[voxel] >= half_max) {
        left_above_half.y = y;
        break;
      }
    }
    for (int y = img.height - 1; y >= max_voxel.y; --y) {
      Voxel voxel(max_voxel.x, y, max_voxel.z);
      if (img[voxel] >= half_max) {
        right_above_half.y = y;
        break;
      }
    }
  }

#if _OPENMP
#pragma omp task shared(img, left_above_half, right_above_half)
#endif
  {
    for (int z = 0; z <= max_voxel.z; ++z) {
      Voxel voxel(max_voxel.x, max_voxel.y, z);
      if (img[voxel] >= half_max) {
        left_above_half.z = z;
        break;
      }
    }
    for (int z = img.depth - 1; z >= max_voxel.z; --z) {
      Voxel voxel(max_voxel.x, max_voxel.y, z);
      if (img[voxel] >= half_max) {
        right_above_half.z = z;
        break;
      }
    }
  }

#if _OPENMP
#pragma omp taskwait
#endif
}

inline void calculate_psf(const VoxelMap& img,
                          const Voxel max_voxel,
                          const F max,
                          const Voxel left_above_half,
                          const Voxel right_above_half,
                          Vector& left,
                          Vector& right,
                          Vector& psf) {
  const auto half_max = max / 2;

  if (left_above_half.x == 0 || right_above_half.x == img.width - 1) {
    psf.x = 0;
    left.x = left_above_half.x;
    right.x = right_above_half.x;
  } else {
    Voxel la(left_above_half.x - 0, max_voxel.y, max_voxel.z);
    Voxel lb(left_above_half.x - 1, max_voxel.y, max_voxel.z);
    Voxel ra(right_above_half.x + 0, max_voxel.y, max_voxel.z);
    Voxel rb(right_above_half.x + 1, max_voxel.y, max_voxel.z);
    left.x = left_above_half.x - (img[la] - half_max) / (img[la] - img[lb]);
    right.x = right_above_half.x + (img[ra] - half_max) / (img[ra] - img[rb]);
    psf.x = right.x - left.x;
  }

  if (left_above_half.y == 0 || right_above_half.y == img.height - 1) {
    psf.y = 0;
    left.y = left_above_half.y;
    right.y = right_above_half.y;
  } else {
    Voxel la(max_voxel.x, left_above_half.y - 0, max_voxel.z);
    Voxel lb(max_voxel.x, left_above_half.y - 1, max_voxel.z);
    Voxel ra(max_voxel.x, right_above_half.y + 0, max_voxel.z);
    Voxel rb(max_voxel.x, right_above_half.y + 1, max_voxel.z);
    left.y = left_above_half.y - (img[la] - half_max) / (img[la] - img[lb]);
    right.y = right_above_half.y + (img[ra] - half_max) / (img[ra] - img[rb]);
    psf.y = right.y - left.y;
  }

  if (left_above_half.z == 0 || right_above_half.z == img.depth - 1) {
    psf.z = 0;
    left.z = left_above_half.z;
    right.z = right_above_half.z;
  } else {
    Voxel la(max_voxel.x, max_voxel.y, left_above_half.z - 0);
    Voxel lb(max_voxel.x, max_voxel.y, left_above_half.z - 1);
    Voxel ra(max_voxel.x, max_voxel.y, right_above_half.z + 0);
    Voxel rb(max_voxel.x, max_voxel.y, right_above_half.z + 1);
    left.z = left_above_half.z - (img[lb] - half_max) / (img[la] - img[lb]);
    right.z = right_above_half.z + (img[rb] - half_max) / (img[ra] - img[rb]);
    psf.z = right.z - left.z;
  }
}

void print_psf(const std::string& fn,
               const VoxelMap& img,
               const F s_pixel,
               std::ostream& out) {
  (void)img;
  Voxel max_voxel;
  F max;
  find_max(img, max_voxel, max);
  Voxel left_above_half, right_above_half;
  find_left_right_above_half(
      img, max_voxel, max, left_above_half, right_above_half);
  Vector left, right, psf;
  calculate_psf(
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
