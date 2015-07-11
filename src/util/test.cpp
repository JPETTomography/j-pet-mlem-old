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
///  - 2d/barrel/event/set
///  - 2d/barrel/detector_set/math
///  - 2d/barrel/detector_set/detect
///  - 2d/barrel/lor/ctor
///  - 2d/barrel/lor/iterator
///  - 2d/barrel/pix_major_system_matrix/ctor
///  - 2d/barrel/pix_major_system_matrix/add
///  - 2d/barrel/pix_major_system_matrix/add_twice
///  - 2d/barrel/pix_major_system_matrix/add_to_all
///  - 2d/barrel/pix_major_system_matrix/to_sparse
///  - 2d/barrel/scanner/math
///  - 2d/barrel/scanner_builder/single_ring/symmetry
///  - 2d/barrel/scanner_builder/multi_ring/symmetry
///  - 2d/barrel/sparse_matrix/symmetric_lor
///  - 2d/barrel/square_detector/intersection
///  - 2d/geometry/circle/init
///  - 2d/geometry/circle/secant
///  - 2d/geometry/circle/secant/math
///  - 2d/geometry/ellipse
///  - 2d/geometry/event
///  - 2d/geometry/line_drawing
///  - 2d/geometry/line_segment
///  - 2d/geometry/pixel_grid/pixel_at
///  - 2d/geometry/pixel_grid/point_at
///  - 2d/geometry/polygon/center
///  - 2d/geometry/polygon/intersection
///  - 2d/geometry/polygon/intersection/math
///  - 2d/geometry/rectangle
///  - 2d/geometry/rectangle/point_generator
///  - 2d/strip/event/conversions1
///  - 2d/strip/event/conversions2
///  - 2d/strip/sensitivity/square
///  - 2d/strip/sensitivity/non_square
///  - 2d/strip/kernel/ctor2
///  - 2d/strip/kernel/bbox
///  - 2d/strip/scanner/pixel_at
///  - 2d/strip/scanner/pixel_center
///  - 3d/geometry/event_generator/spherical_distribution
///  - 3d/geometry/event_generator/voxel_event_generator
///  - 3d/geometry/event_generator/cylinder_event_generator
///  - 3d/geometry/event_generator/ball_event_generator
///  - 3d/geometry/event_generator/ellipsoid_event_generator
///  - 3d/geometry/matrix/initialisation
///  - 3d/geometry/matrix/vector_multiplication
///  - 3d/geometry/matrix/arithmetic_assignment_operators
///  - 3d/geometry/matrix/arithmetic_operators
///  - 3d/geometry/matrix/scalar_multiplication
///  - 3d/geometry/phantom_builder/rapid_json
///  - 3d/geometry/phantom_builder/angular_distribution
///  - 3d/geometry/phantom_builder/angular_distribution/spherical
///  - 3d/geometry/phantom_builder/phantom
///  - 3d/geometry/phantom/cylinder_region
///  - 3d/geometry/phantom/cylinder
///  - 3d/geometry/phantom/ellipsoid
///  - 3d/geometry/point/init
///  - 3d/geometry/point/arithmetic assignemt
///  - 3d/geometry/point/difference
///  - 3d/geometry/point/distance
///  - 3d/geometry/point/nearest_distance
///  - 3d/geometry/vector/init
///  - 3d/geometry/vector/arithmetic_assignement
///  - 3d/geometry/vector/arithmetics
///  - 3d/geometry/vector/logical
///  - 3d/geometry/voxel_grid
///  - 3d/geometry/voxel_set
///  - 3d/hybrid/detector_set/escape_through_endcap
///  - 3d/hybrid/detector_set/detect
///  - 3d/hybrid/sensitivity_mapper
///  - common/phantom_monte_carlo/point_source
///  - common/phantom_monte_carlo/phantom_region
///  - util/array
///  - util/grapher/detector
///  - util/grapher/big_barrel
///  - util/grapher/big_barrel/lor
///  - util/grapher/big_barrel/segment
///  - util/grapher/big_barrel/circle
///  - util/grapher/big_barrel/pixel
///  - util/random
///  - util/sort
/// 82 test cases

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_USE_ANSI_COLOUR_CODES
#include "catch.hpp"
