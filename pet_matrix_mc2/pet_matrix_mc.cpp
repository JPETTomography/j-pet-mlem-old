// PET System Matrix Calculator
// Author: Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//
// Using Monte Carlo method and square scintilators.

#include <random>
#include <iostream>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
namespace po = boost::program_options;

using namespace std;

template <typename F> F deg(F rad) { return rad * 180/M_PI; }
template <typename F> F rad(F deg) { return deg * M_PI/180; }

typedef map<pair<int, int>, vector<int>> lor_map;

// produces secant points using circle/line intersection as a equation system solution
// degenerated  http://www.wolframalpha.com/input/?i=solve+x%5E2%2By%5E2%3DR%5E2%2C+x%3DC+for+x%2C+y
// normal slope http://www.wolframalpha.com/input/?i=solve+x%5E2%2By%5E2%3DR%5E2%2C+y%3DM*x%2BB+for+x%2C+y
template <typename F = double>
class scintilator_ring {
public:
  scintilator_ring(F radious)
  : r(radious)
  , r2(radious*radious) {}

  typedef F angle_type;
  typedef pair<F, F> point_type;
  typedef pair<point_type, point_type> secant_type;

  secant_type secant(point_type &p, angle_type phi) {
    // switch coordinates when close to M_PI_2 to get accurate results
    bool degen = phi > M_PI_4 && phi < 3.0 * M_PI_4;
    auto m = !degen ? tan(phi) : tan(phi - M_PI_2);
    auto b = !degen ? (p.first - m * p.second) : (p.second - m * p.first);
    auto m2 = m*m;
    auto b2 = m*m;
    auto sq = sqrt(-b2 + m2*r2 + r2);

    auto s = make_pair(
      make_pair((-sq - b*m ) / (m2 + 1.0), (b - m*sq) / (m2 + 1.0)),
      make_pair(( sq - b*b ) / (m2 + 1.0), (b + m*sq) / (m2 + 1.0))
    );
    if (degen) return make_pair(
        make_pair(s.first.second,  s.first.first),
        make_pair(s.second.second, s.second.first)
      );
    return s;
  }
private:
  F r;
  F r2;
};

int main(int argc, char *argv[]) {

  size_t n_pixels;
  size_t n_detectors;
  size_t n_emissions;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "shows this help")
    ("n-pixels,n",    po::value<size_t>(&n_pixels)->default_value(256),
     "set number of pixels in one dimension")
    ("n-detectors,d", po::value<size_t>(&n_detectors)->default_value(64),
     "set number of ring detector")
    ("n-emissions,e", po::value<size_t>(&n_emissions)->default_value(1),
     "emissions per pixel")
  ;

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch(const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  lor_map mc;
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> ydis(0, n_pixels/2-1);
  uniform_real_distribution<> rdis(0, 1);
  uniform_real_distribution<> hdis(0, .5);
  uniform_real_distribution<> phidis(0, M_PI);

  auto radious = sqrt(n_pixels * n_pixels / 2);
  cout << "ring radious = " << radious << endl;
  scintilator_ring<> sr(radious);

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

#if 1
    scintilator_ring<>::point_type p(
      px + rdis(gen),
      py + ((px == py) ? hdis(gen) : rdis(gen))
    );
    auto phi = phidis(gen);
#else
    scintilator_ring<>::point_type p(0.0, 0.0);
    auto phi = M_PI/2;
#endif
    scintilator_ring<>::secant_type s = sr.secant(p, phi);

    cout << '(' << p.first << ',' << p.second << ')'
           << ' ' << deg(phi)
         << " s1 = (" << s.first.first  << ',' << s.first.second  << ')'
           << ' ' << deg(atan2(s.first.second, s.first.first))
         << " s2 = (" << s.second.first << ',' << s.second.second << ')'
           << ' ' << deg(atan2(s.second.second, s.second.first))
         << endl;
  }

  return 0;
}
