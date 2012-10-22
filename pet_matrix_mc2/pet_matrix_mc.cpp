// PET System Matrix Calculator
// Author: Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//
// Using Monte Carlo method and square scintilators.

#include <random>
#include <iostream>
#include <cmdline.h>

using namespace std;

#include "circle.h"

template <typename F> F deg(F rad) { return rad * 180/M_PI; }
template <typename F> F rad(F deg) { return deg * M_PI/180; }

typedef map<pair<int, int>, vector<int>> lor_map;

int main(int argc, char *argv[]) {

  cmdline::parser cl;

  cl.add<size_t>("n-pixels",    'n', "number of pixels in one dimension", false, 256);
  cl.add<size_t>("n-detectors", 'd', "number of ring detectors",          false, 64);
  cl.add<size_t>("n-emissions", 'e', "emissions per pixel",               false, 1);

  cl.parse_check(argc, argv);

  size_t n_pixels    = cl.get<size_t>("n-pixels");
  size_t n_detectors = cl.get<size_t>("n-detectors");
  size_t n_emissions = cl.get<size_t>("n-emissions");

  lor_map mc;
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> ydis(0, n_pixels/2-1);
  uniform_real_distribution<> rdis(0, 1);
  uniform_real_distribution<> hdis(0, .5);
  uniform_real_distribution<> phidis(0, M_PI);

  auto radious = sqrt(n_pixels * n_pixels / 2);
  cout << "circle radious = " << radious << endl;
  circle<> c(radious);

  for (auto n = 0; n < n_emissions; ++n) {
    // generate random y pixel
    auto py = ydis(gen);
    // generate random x pixel <= y
    uniform_int_distribution<> xdis(0, py);
    auto px = xdis(gen);

    auto lor = lor_map::key_type(0, n_detectors/2);
    auto lor_pixels = mc[lor];

    lor_pixels.reserve(n_pixels * (n_pixels+1)/2);
    lor_pixels[py * (py+1)/2 + px] ++;

#if 0
    decltype(c)::point_type p(
      px + rdis(gen),
      py + ((px == py) ? hdis(gen) : rdis(gen))
    );
    auto phi = phidis(gen);
#else
    decltype(c)::point_type p(0.0, 0.0);
    auto phi = M_PI/2;
#endif
    auto s = c.secant(p, phi);

    cout << '(' << p.x << ',' << p.y << ')'
           << ' ' << deg(phi)
         << " s1 = (" << s.first.x  << ',' << s.first.y  << ')'
           << ' ' << deg(atan2(s.first.y, s.first.x))
         << " s2 = (" << s.second.x << ',' << s.second.y << ')'
           << ' ' << deg(atan2(s.second.y, s.second.x))
         << endl;
  }

  return 0;
}
