#pragma once

#include <cstring>

#include "../geometry/voxel.h"
#include "../geometry/voxel_map.h"
#include "../geometry/vector.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

namespace PET3D {
/// Generic tools
namespace Tool {

/// Generates PSF FWHM out of VoxelMap
struct PSF {
  using Voxel = PET3D::Voxel<S>;
  using VoxelMap = PET3D::VoxelMap<Voxel, F>;
  using Vector = PET3D::Vector<F>;

  static void find_max(const VoxelMap& img,
                       Voxel& max_voxel,
                       F& max,
                       S padding = 0) {
    auto thread_max_voxels = new (alloca(sizeof(Voxel) * omp_get_max_threads()))
        Voxel[omp_get_max_threads()];
    auto thread_maxes = new (alloca(sizeof(F) * omp_get_max_threads()))
        F[omp_get_max_threads()]();
    auto x_padding = std::min(padding, S(img.width - 1));
    auto y_padding = std::min(padding, S(img.height - 1));
    auto z_padding = std::min(padding, S(img.depth - 1));
#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (S z = z_padding; z < img.depth - z_padding; ++z) {
      Voxel l_max_voxel = thread_max_voxels[omp_get_thread_num()];
      auto l_max = thread_maxes[omp_get_thread_num()];
      for (S y = y_padding; y < img.height - y_padding; ++y) {
        for (S x = x_padding; x < img.width - x_padding; ++x) {
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
      if (thread_maxes[t] > max ||
          (thread_maxes[t] == max && thread_max_voxels[t] < max_voxel)) {
        max = thread_maxes[t];
        max_voxel = thread_max_voxels[t];
      }
    }
  }

  static void find_left_right_above_half(const VoxelMap& img,
                                         const Voxel max_voxel,
                                         const F max,
                                         Voxel& left_above_half,
                                         Voxel& right_above_half) {
    const auto half_max = max / 2;
#if _OPENMP && !_MSC_VER
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

#if _OPENMP && !_MSC_VER
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

    if (img.depth < 2) {
      left_above_half.z = 0;
      right_above_half.z = 0;
#if _OPENMP && !_MSC_VER
#pragma omp taskwait
#endif
      return;
    }

#if _OPENMP && !_MSC_VER
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

#if _OPENMP && !_MSC_VER
#pragma omp taskwait
#endif
  }

  static void calculate(const VoxelMap& img,
                        const Voxel max_voxel,
                        const F max,
                        const Voxel left_above_half,
                        const Voxel right_above_half,
                        Vector& left,
                        Vector& right,
                        Vector& psf) {
    const auto half_max = max / 2;

    if (left_above_half.x == 0 || right_above_half.x == img.width - 1) {
      psf.x = -1;
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
      psf.y = -1;
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

    if (img.depth < 2) {
      psf.z = -1;
      left.z = 0;
      right.z = 0;
      return;
    }

    if (left_above_half.z == 0 || right_above_half.z == img.depth - 1) {
      psf.z = -1;
      left.z = left_above_half.z;
      right.z = right_above_half.z;
    } else {
      Voxel la(max_voxel.x, max_voxel.y, left_above_half.z - 0);
      Voxel lb(max_voxel.x, max_voxel.y, left_above_half.z - 1);
      Voxel ra(max_voxel.x, max_voxel.y, right_above_half.z + 0);
      Voxel rb(max_voxel.x, max_voxel.y, right_above_half.z + 1);
      left.z = left_above_half.z - (img[la] - half_max) / (img[la] - img[lb]);
      right.z = right_above_half.z + (img[ra] - half_max) / (img[ra] - img[rb]);
      psf.z = right.z - left.z;
    }
  }
};

}  // Hybrid
}  // PET3D
