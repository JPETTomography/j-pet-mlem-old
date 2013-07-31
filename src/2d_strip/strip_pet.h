#ifndef STRIP_PET_H
#define STRIP_PET_H

#include <cmath>
#include <vector>
#include <fstream>
#include <ctime>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/additive_combine.hpp>
#include <boost/random/normal_distribution.hpp>



typedef boost::ecuyer1988 base_generator_type;
typedef boost::normal_distribution<> distribution_type;
typedef boost::variate_generator<base_generator_type&, distribution_type>
gen_type;

template <typename T = float> class Strip_PET {
 private:
  T R_distance;
  T Scentilator_length;
  T _a, _b, _x, _y, _phi;
  T _sin;
  T _cos;
  T _inv_a2;
  T _inv_b2;
  std::vector<event<T>> event_list;
  std::vector<scintillator<>> scientilator_list;

 public:
  Strip_PET(T& R_distance,
            T& Scentilator_length,
            T& x,
            T& y,
            T& a,
            T& b,
            T& phi)
      : R_distance(R_distance),
        Scentilator_length(Scentilator_length),
        _x(x),
        _y(y),
        _a(a),
        _b(b),
        _phi(phi) {
    _sin = sin(static_cast<double>(_phi));
    _cos = cos(static_cast<double>(_phi));
    _inv_a2 = 1.0 / (_a * _a);
    _inv_b2 = 1.0 / (_b * _b);

    scientilator_list.push_back(
        scintillator<>(R_distance, 0.0f, Scentilator_length));
    scientilator_list.push_back(
        scintillator<>(-R_distance, 0.0f, Scentilator_length));
  }

  bool in(T x, T y) const {

    T dx = x - _x;
    T dy = y - _y;

    T r_x = dx * _cos + dy * _sin;
    T r_y = -dx * _sin + dy * _cos;

    T r2 = r_x * r_x * _inv_a2 + r_y * r_y * _inv_b2;

    return r2 < 1.0;
  }

  void emit_event() {

    event<T> temp_event;

    base_generator_type generator(42);

    gen_type normal_distribution(generator, distribution_type(0.0, 10.0f));

    T rand_y = 2.0;
    T rand_z = 2.0;
    T angle = (M_PI / 180);

    for (int i = 0; i < 90; ++i) {
      if (in(rand_y, rand_z)) {

        angle = (3.14 / 180) * i;

        T z_u =
            rand_z + (R_distance - rand_y) * tan(angle) + normal_distribution();
        T z_d =
            rand_z - (R_distance + rand_y) * tan(angle) + normal_distribution();
        T dl = -static_cast<T>(2) *
                   sqrt(static_cast<T>(1) + (tan(angle) * tan(angle))) +
               normal_distribution();

        if (std::abs(z_u) < scientilator_list[0].get_l() &&
            std::abs(z_d) < scientilator_list[1].get_l()) {

          temp_event.z_u = z_u;
          temp_event.z_d = z_d;
          temp_event.dl = dl;

          std::cout << "EVENT: " << angle << " " << z_u << " " << z_d << " "
                    << dl << std::endl;
          event_list.push_back(temp_event);
        }
      }
    }
  }
};

#endif  // STRIP_PET_H
