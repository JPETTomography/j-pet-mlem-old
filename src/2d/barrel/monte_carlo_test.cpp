#include "util/test.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/geometry/phantom.h"
#include "2d/geometry/pixel_map.h"
#include "common/phantom_monte_carlo.h"

#include "common/model.h"
#include "common/types.h"

using RNG = std::mt19937;
using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner = PET2D::Barrel::GenericScanner<Detector, short, 64>;
using Phantom = PET2D::Phantom<RNG, F>;
using Allways = Common::AlwaysAccept<F>;
using Scintillator = Common::ScintillatorAccept<F>;
using Pixel = PET2D::Pixel<S>;
using Point = PET2D::Point<F>;
using Event = PET2D::Event<F>;
using Image = PET2D::PixelMap<Pixel, Hit>;

using MonteCarlo = Common::PhantomMonteCarlo<Phantom, Scanner, Image>;

namespace {
F strip_width = 0.005;
F strip_height = 0.019;
F strip_distance = 0.410;
F inner_radius = (strip_distance - strip_height) / 2;
F strip_length = 0.300;
}

TEST("2d/barrel/phantom_monte_carlo") {
    auto scanner2d = PET2D::Barrel::ScannerBuilder<Scanner>::build_single_ring(
        inner_radius, 32, strip_width, strip_height);
}
