#pragma once

#include <vector>
#include <sstream>
#include <iomanip>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/progress.h"
#include "util/delegate.h"
#include "common/types.h"

#include "../../geometry/point.h"
#include "../../geometry/pixel.h"
#include "../../geometry/pixel_map.h"

#include "../response.h"
#include "../scanner.h"

namespace PET2D {
namespace Strip {
/// CUDA optimized subimplementation
namespace GPU {
namespace Reconstruction {

using Point = PET2D::Point<F>;
using Pixel = PET2D::Pixel<S>;
using Output = PixelMap<Pixel, F>;

/// CUDA entry-point function
template <typename F>
void run(Scanner<F, S>& scanner,
         Response<F>* responses,
         int n_responses,
         Output& input_output,
         int n_iteration_blocks,
         int n_iterations_in_block,
         util::delegate<void(int iteration, const Output& output)> output,
         util::delegate<void(int completed, bool finished)> progress,
         int device,
         int n_blocks,
         int n_threads_per_block,
         util::delegate<void(const char* device_name)> device_name);

}  // Reconstruction
}  // GPU
}  // Strip
}  // PET2D
