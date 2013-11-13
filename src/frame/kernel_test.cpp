#include "catch.hpp"

#include "detector_frame.h"
#include "error_model.h"
#include "kernel.h"
#include "2d_xy/square_detector.h"

typedef float FLOAT;

const FLOAT degree = M_PI / 180.0;

const FLOAT R = 450.0;
const FLOAT L = 300.0;
const FLOAT sigma_z = 10.0;
const FLOAT sigma_l = 63.0;

typedef Kernel<FLOAT, SquareDetector, ConstantErrorsModel> KernelType;
typedef KernelType::EventKernelType EventKernelType;

#if TEST_FRAME_KERNEL
// FIXME: this doesn't work as it is NYI, disable
TEST_CASE("Kernel/Create", "Create") {

  KernelType kernel(Detector<FLOAT>(R, L),
                    ConstantErrorsModel<FLOAT>(sigma_z, sigma_l));

  EventImageAngle<FLOAT> event(300.0, 0.0, 18.0 * Degree);

  EventKernelType event_kernel = kernel.MakeEventKernel(event);

  REQUIRE(std::abs(event_kernel(0.0, 0.0) - 9.5026156e-8) < 1e-12);
}
#endif
