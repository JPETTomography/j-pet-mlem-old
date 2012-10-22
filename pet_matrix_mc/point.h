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

  // I know it is bad idea to count all over again
  // sin/cos for given point, but this will be used
  // only for initialization.
  point rotated(F phi) {
    auto sin_phi = std::sin(phi);
    auto cos_phi = std::cos(phi);
    return {x*cos_phi - y*sin_phi, x*sin_phi + y*cos_phi};
  }

  point operator + (point &p) {
    return {x+p.x, y+p.y};
  }
};
