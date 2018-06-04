#include "util/test.h"

#include <cmath>

#include "geometry.h"

#include "common/types.h"

using Geometry = PET2D::Barrel::Geometry<F, S>;
using Grid = PET2D::PixelGrid<F, S>;
using Point = PET2D::Point<F>;

TEST("2d/geometry/geometry") {}

TEST("2d/geometry/geometry/write_read") {
  Grid grid(192, 128, 0.1, Point(-1, -1));
  Geometry geometry(48, grid);
  {
    auto fn = "geometry_test.txt"_temp;
    std::ofstream out(fn);
    out << geometry;
    out.close();

    std::ifstream in(fn);
    REQUIRE(in);
    Geometry read_geometry(in);
    in.close();

    CHECK(geometry.n_detectors == read_geometry.n_detectors);
    CHECK(geometry.grid.n_columns == read_geometry.grid.n_columns);
    CHECK(geometry.grid.n_rows == read_geometry.grid.n_rows);
    CHECK(geometry.grid.pixel_size == read_geometry.grid.pixel_size);

    REQUIRE(std::remove(fn.c_str()) == 0);
  }
  {
    auto fn = "geometry_test.bin"_temp;
    util::obstream out(fn);
    out << geometry;
    out.close();

    util::ibstream in(fn);
    REQUIRE(in);
    Geometry read_geometry(in);
    in.close();

    CHECK(geometry.n_detectors == read_geometry.n_detectors);
    CHECK(geometry.grid.n_columns == read_geometry.grid.n_columns);
    CHECK(geometry.grid.n_rows == read_geometry.grid.n_rows);
    CHECK(geometry.grid.pixel_size == read_geometry.grid.pixel_size);

    REQUIRE(std::remove(fn.c_str()) == 0);
  }
}
