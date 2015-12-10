/// \page cmd_3d_tool_crop 3d_tool_crop
/// \brief Crops image into subimage
///
/// Usage
/// -----
/// \verboutput 3d_tool_crop
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

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  PET3D::Tool::add_crop_options(cl);
  cl.parse_check(argc, argv);
  PET3D::Tool::calculate_psf_options(cl, argc);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  auto n_pixels = cl.get<int>("n-pixels");
  auto s_pixel = cl.get<double>("s-pixel");
  PixelGrid pixel_grid(n_pixels, n_pixels, s_pixel);

  auto output_name = cl.get<cmdline::path>("output");
  auto output_base_name = output_name.wo_ext();
  auto output_ext = output_name.ext();
  auto output_txt = output_ext == ".txt";

  int n_planes = cl.get<int>("n-planes");
  VoxelGrid grid(pixel_grid, -s_pixel * n_planes / 2, n_planes);

  std::cerr << "   voxel grid = "  // grid size:
            << grid.pixel_grid.n_columns << " x " << grid.pixel_grid.n_rows
            << " x " << grid.n_planes << std::endl;
  std::cerr << "   voxel size = " << s_pixel << std::endl;

  auto crop_size = cl.get<int>("crop");
  auto crop_origin = PET3D::Voxel<S>(
      cl.get<int>("crop-x"), cl.get<int>("crop-y"), cl.get<int>("crop-z"));
  VoxelMap cropped_img(crop_size, crop_size, crop_size);

  std::cerr << "  crop origin = "
            << "(" << crop_origin.x << ", " << crop_origin.y << ", "
            << crop_origin.z << ")" << std::endl;

  for (const auto& in_fn : cl.rest()) {
    if (cmdline::path(in_fn).ext() == ".txt") {
      throw("text files are not supported by this tool: " + in_fn);
    }
    std::cerr << "        image = " << in_fn << std::endl;
#if USE_MMAP
    auto fd = open(in_fn.c_str(), O_RDONLY);
    if (fd == -1) {
      throw("cannot open: " + in_fn);
    }
    const auto data_size = grid.n_voxels * sizeof(F);
    F* data = (F*)mmap(NULL, data_size, PROT_READ, MAP_SHARED, fd, 0);
    VoxelMap img(
        grid.pixel_grid.n_columns, grid.pixel_grid.n_rows, grid.n_planes, data);
#else
    VoxelMap img(
        grid.pixel_grid.n_columns, grid.pixel_grid.n_rows, grid.n_planes);
    util::ibstream bin(in_fn);
    if (!bin.is_open()) {
      throw("cannot open: " + fn);
    }
    bin >> img;
#endif

    auto index = cmdline::path(in_fn).scan_index();
    auto fn = index >= 0 ? output_base_name.add_index(index, index)
                         : output_base_name;

    std::cerr << "       output = " << fn + output_ext << std::endl;
    cropped_img.copy(img, crop_origin);
    util::nrrd_writer nrrd(fn + ".nrrd", fn + output_ext, output_txt);
    nrrd << cropped_img;
    if (output_txt) {
      std::ofstream txt(fn + ".txt");
      txt << cropped_img;
    } else {
      util::obstream bin(fn + output_ext);
      bin << cropped_img;
    }

#if USE_MMAP
    munmap(data, data_size);
    close(fd);
#endif
  }

  CMDLINE_CATCH
}
