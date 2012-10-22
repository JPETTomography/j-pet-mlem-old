#pragma once

template <typename F = double>
struct point {
  point(F _x, F _y)
  : x(_x)
  , y(_y) {}

  point(std::pair<F, F> p)
  : x(p.first)
  , y(p.second) {}

  F x, y;
};
