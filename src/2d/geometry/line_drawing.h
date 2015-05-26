#ifndef LINE_DRAWING
#define LINE_DRAWING

#include <limits>

#include "2d/geometry/pixel_grid.h"
namespace PET2D {

/*
 *  A Fast Voxel Traversal Algorithm for Ray Tracing
 *  John Amanatides
 *  Andrew Woo
 */

template <typename Grid, typename I>
void draw_line(const Point<typename Grid::F>& start,
               const Point<typename Grid::F>& end,
               const Grid grid,
               I pixels) {

  using Vector = typename Grid::Vector;

  using Pixel = typename Grid::Pixel;
  using Point = typename Grid::Point;
  using F = typename Grid::F;
  using S = typename Grid::S;

  typename Grid::Vector diff = end - start;

  F length = diff.length();
  F s = diff.y / length;
  F c = diff.x / length;

  Pixel start_pix = grid.pixel_at(start);
  Pixel end_pix = grid.pixel_at(end);

  typename Grid::S step_x = (end.x >= start.x) ? 1 : -1;
  typename Grid::S step_y = (end.y >= start.y) ? 1 : -1;

  Point start_pixel_center = grid.center_at(start_pix.x, start_pix.y);
  F x_border = start_pixel_center.x + step_x * grid.pixel_size / 2;
  F y_border = start_pixel_center.y + step_y * grid.pixel_size / 2;

  F t_max_x;
  if (start.x != end.x) {
    t_max_x = (x_border - start.x) / c;
  } else {
    t_max_x=step_x*std::numeric_limits<F>::max();
  }

  F t_max_y;
  if (start.y != end.y) {
    t_max_y = (y_border - start.y) / s;
  } else {
    t_max_y = step_y*std::numeric_limits<F>::max();
  }

  F t_delta_x = grid.pixel_size / c;
  F t_delta_y = grid.pixel_size / s;


  S ix=start_pix.x, iy=start_pix.y;
  pixels=Pixel(ix,iy);
  pixels++;
  while (ix != end_pix.x || iy!=end_pix.y) {
    if (t_max_x < t_max_y) {
      t_max_x+=t_delta_x;
      ix+=step_x;
    } else {
      t_max_y+=t_delta_y;
      iy+=step_y;
    }
    pixels=Pixel(ix,iy);
    pixels++;
  }
}
}
#endif  // LINE_DRAWING
