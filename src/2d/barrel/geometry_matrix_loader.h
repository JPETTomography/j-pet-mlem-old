#ifndef GEOMETRY_MATRIX_LOADER
#define GEOMETRY_MATRIX_LOADER

#include <string>

#include "2d/barrel/lor_geometry.h"
#include "2d/barrel/geometry.h"
#include "2d/barrel/sparse_matrix.h"

template <typename F, typename S, typename Hit>
void load_system_matrix_from_file(std::string system_matrix_file_name,
                                  PET2D::Barrel::Geometry<F, S>& geometry,
                                  bool verbose) {
  using Point = PET2D::Point<F>;
  using Pixel = PET2D::Pixel<S>;
  using LOR = PET2D::Barrel::LOR<S>;
  // geometry.erase_pixel_info();

  util::ibstream in_matrix(system_matrix_file_name);
  if (!in_matrix.is_open()) {
    throw("cannot open system matrix file: " + system_matrix_file_name);
  }
  PET2D::Barrel::SparseMatrix<Pixel, LOR, Hit> matrix(in_matrix);
  if (verbose) {
    std::cout << "read in system matrix: " << system_matrix_file_name
              << std::endl;
  }
  matrix.sort_by_lor_n_pixel();
  geometry.sort_all_by_index();
  // matrix.merge_duplicates();
  F n_emissions = F(matrix.n_emissions());
  if (geometry.grid.n_columns != matrix.n_pixels_in_row()) {
    throw("mismatch in number of pixels with matrix");
  }
  if (matrix.triangular()) {
    throw("matrix is not full");
  }

  for (auto& element : matrix) {
    auto lor = element.lor;
    F weight = element.hits / n_emissions;
    geometry.put_pixel(lor, element.pixel, weight);
  }
  geometry.sort_all();
}

#endif  // GEOMETRY_MATRIX_LOADER
