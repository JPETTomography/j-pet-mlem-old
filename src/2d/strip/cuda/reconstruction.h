#pragma once

#include <vector>
#include <sstream>
#include <iomanip>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/progress.h"
#include "util/delegate.h"

#include "../response.h"
#include "../scanner.h"

namespace PET2D {
namespace Strip {
/// CUDA optimized subimplementation
namespace GPU {
namespace Reconstruction {

/// CUDA entry-point function
template <typename F>
void run(Scanner<F, short>& scanner,
         Response<F>* responses,
         int n_responses,
         int n_iteration_blocks,
         int n_iterations_in_block,
         util::delegate<void(int, F*)> output,
         util::delegate<void(int, bool)> progress,
         int device,
         int n_blocks,
         int n_threads_per_block,
         util::delegate<void(const char*)> device_name);

}  // Reconstruction
}  // GPU
}  // Strip
}  // PET2D
