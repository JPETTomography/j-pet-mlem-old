/// \page test Testing
/// \brief More about unit-testing in this project
///
/// [catch]: https://github.com/philsquared/Catch
///
/// This project uses [Catch][catch] for unit-testing.
///
/// Tests are not build by default. In order to build and run tests run:
///
///     make test && ./test
///
/// Files containing tests have `_test.cpp` suffix.
///
/// Usage
/// -----
/// \verbinclude src/util/test.txt
///
/// All available test cases
/// ------------------------
///  - 2d/barrel/circle_detector/ctor
///  - 2d/barrel/circle_detector/move
///  - 2d/barrel/circle_detector/intersection
///  - 2d/barrel/scanner/math
///  - 2d/barrel/detector_set/math
///  - 2d/barrel/detector_set/detect
///  - 2d/barrel/event/set
///  - 2d/barrel/lor/ctor
///  - 2d/barrel/lor/iterator
///  - 2d/barrel/pix_major_system_matrix/ctor
///  - 2d/barrel/pix_major_system_matrix/add
///  - 2d/barrel/pix_major_system_matrix/add_twice
///  - 2d/barrel/pix_major_system_matrix/add_to_all
///  - 2d/barrel/pix_major_system_matrix/to_sparse
///  - 2d/barrel/phantom/phantom
///  - 2d/barrel/square_detector/intersection
///  - 2d/geometry/circle/init
///  - 2d/geometry/circle/secant
///  - 2d/geometry/circle/secant/math
///  - 2d/geometry/ellipse/elliptical_region
///  - 2d/geometry/polygon/center
///  - 2d/geometry/polygon/intersection
///  - 2d/geometry/polygon/intersection/math
///  - 2d/strip/detector/pixel_at
///  - 2d/strip/detector/pixel_center
///  - 2d/strip/event/conversions1
///  - 2d/strip/event/conversions2
///  - 2d/strip/sensitivity/square
///  - 2d/strip/sensitivity/non_square
///  - 2d/strip/kernel/ctor2
///  - 2d/strip/kernel/bbox
///  - geometry/point/2d
///  - util/array
///  - util/random
///
/// 34 test cases

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_USE_ANSI_COLOUR_CODES
#include "catch.hpp"
