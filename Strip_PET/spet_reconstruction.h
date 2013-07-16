#ifndef SPET_RECONSTRUCTION_
#define SPET_RECONSTRUCTION_H

#include <cmath>
#include <vector>
#include <algorithm>
#include "data_structures.h"

// algebra
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

#define INVERSE_PI 0.31830988
#define INVERSE_TWO_PI 0.15915494

using namespace boost::numeric::ublas;

template <typename T = float>
class spet_reconstruction {

  typedef std::pair<int, int> pixel_location;

 private:
  T R_distance;
  T Scentilator_length;
  int n_pixels;
  int n_pixels_2;
  T pixel_size;
  T sigma_z;
  T sigma_dl;
  T gamma;
  std::vector<std::vector<event<T> > > event_list;
  std::vector<std::vector<T> > inverse_correlation_matrix;
  T inverse_correlation_matrix_factor;
  matrix<T> C;
  matrix<T> inverse_C;
  std::vector<T> rho;
  // vector elements
  T two;
  T one;
  T zero;
  T half;

 public:
  spet_reconstruction(T &R_distance, T &Scentilator_length, int &n_pixels,
                      T &pixel_size, T &sigma_z, T &sigma_dl, T &gamma)
      : R_distance(R_distance),
        Scentilator_length(Scentilator_length),
        n_pixels(n_pixels),
        pixel_size(pixel_size),
        sigma_z(sigma_z),
        sigma_dl(sigma_dl),
        gamma(gamma) {

    n_pixels_2 = n_pixels / 2;

    event_list.resize(n_pixels * n_pixels);
    rho.resize(n_pixels * n_pixels);

    set_inverse_correlation_matrix();

    two = static_cast<T>(2);
    one = static_cast<T>(1);
    zero = static_cast<T>(0);
    half = static_cast<T>(0.5);
  }

  void set_inverse_correlation_matrix() {

    inverse_correlation_matrix.resize(3, std::vector<T>(3, zero));

    T pow_sigma_z = sigma_z * sigma_z;
    T pow_sigma_dl = sigma_dl * sigma_dl;
    T pow_gamma = gamma * gamma;

    inverse_correlation_matrix[0][0] = pow_gamma - (pow_sigma_z * pow_sigma_dl);
    inverse_correlation_matrix[0][1] = pow_gamma;
    inverse_correlation_matrix[0][2] = gamma * pow_sigma_z;

    inverse_correlation_matrix[1][0] = pow_gamma;
    inverse_correlation_matrix[1][1] = pow_gamma - (pow_sigma_z * pow_sigma_dl);
    inverse_correlation_matrix[1][2] = -gamma * pow_sigma_z;

    inverse_correlation_matrix[2][0] = gamma * pow_sigma_z;
    inverse_correlation_matrix[2][1] = -gamma * pow_sigma_z;
    inverse_correlation_matrix[2][2] = -pow_sigma_z * pow_sigma_z;

    inverse_correlation_matrix_factor =
        one / ((pow_gamma * pow_sigma_z) + (pow_gamma * pow_sigma_z) +
               (pow_sigma_z * pow_sigma_z * pow_sigma_dl));
  }

  T multiply_elements(std::vector<T> &vec_a, std::vector<T> &vec_b) {

    T a[vec_a.size()];
    T output = zero;
    // add AVX
    for (unsigned i = 0; i < vec_a.size(); ++i) {
      for (unsigned j = 0; j < vec_a.size(); j++) {

        a[i] += vec_a[j] * inverse_correlation_matrix[j][i];
      }
      output += a[i] * vec_b[i];
    }

    return inverse_correlation_matrix_factor * output;
  }

  void load_input(std::string fn) {

    ibstream in(fn, std::ios::binary);
    event<T> temp_event;

    int number_of_pixels;
    float pixel_s;
    int iter;
    int number_of_event_in_file;

    in >> number_of_pixels;
    in >> pixel_s;
    in >> iter;
    in >> number_of_event_in_file;
#if DEBUG == 1
    std::cout << number_of_pixels << " " << pixel_s << " " << iter << " "
              << number_of_event_in_file << std::endl;
    std::cout << "VECTOR SIZE: " << event_list.size() << std::endl;

#endif
    for (;;) {

      T z_u, z_d, dl;

      in >> z_u >> z_d >> dl;

      temp_event.z_u = z_u;
      temp_event.z_d = z_d;
      temp_event.dl = dl;

      T tan = event_tan(z_u, z_d);
      T y = event_y(z_d, tan);
      T z = event_z(z_u, z_d, y, tan);

      pixel_location p = in_pixel(y, z);

#if DEBUG == 1
      std::cout << "DATA: " << z_u << "  " << z_d << " " << dl << std::endl;
      std::cout << "PIXEL: " << p.first << " " << p.second << " "
                << p.first *n_pixels + p.second << std::endl;
#endif
      event_list[p.first * n_pixels + p.second].push_back(temp_event);

      if (in.eof()) {
        break;
      }
    }

#if DEBUG == 1

    std::cout << "VECTOR" << std::endl;

    typename std::vector<std::vector<event<T> > >::iterator it;
    typename std::vector<event<T> >::iterator jt;

    for (it = event_list.begin(); it != event_list.end(); ++it) {

      for (jt = (*it).begin(); jt != (*it).end(); ++jt) {

        std::cout << (*jt).z_u << " " << (*jt).z_d << " " << (*jt).dl
                  << std::endl;
      }
    }
#endif
  }

  T sensitivity(T &y, T &z) {

    T a = (half * Scentilator_length - z) / (R_distance - y);
    T b = (half * Scentilator_length + z) / (R_distance + y);
    T c = (half * Scentilator_length + z) / (R_distance - y);

    return INVERSE_PI *
           (std::atan(std::min(a, b)) - std::atan(std::min(-c, -b)));
  }

  T kernel(T &y, T &tan, pixel_location &pixel) {

    T angle = std::atan(tan);

    T event_inv_cos = one / std::cos(angle);
    T pow_event_inv_cos = event_inv_cos * event_inv_cos;

    std::vector<T> vec_o(3, zero);
    std::vector<T> vec_a(3, zero);
    std::vector<T> vec_b(3, zero);

    vec_o[0] = -(pixel.first + y - R_distance) * tan * pow_event_inv_cos;
    vec_o[1] = -(pixel.first + y + R_distance) * tan * pow_event_inv_cos;
    vec_o[2] = -(pixel.first + y) * event_inv_cos * (one + (two * tan * tan));

    vec_a[0] = -(pixel.first + y - R_distance) * pow_event_inv_cos;
    vec_a[1] = -(pixel.first + y + R_distance) * pow_event_inv_cos;
    vec_a[2] = -(pixel.first + y) * event_inv_cos * tan;

    vec_b[0] = pixel.second - (pixel.first * tan);
    vec_b[1] = pixel.second - (pixel.first * tan);
    vec_b[2] = -two * pixel.first * event_inv_cos;

    T element_aa = multiply_elements(vec_a, vec_a);
    T element_ba = multiply_elements(vec_b, vec_a);
    T element_bb = multiply_elements(vec_b, vec_b);
    T element_ob = multiply_elements(vec_o, vec_b);
    T element_sum_aa_ob = element_aa + (two * element_ob);

    T element_before_exp = INVERSE_TWO_PI * (one / sqrt(element_sum_aa_ob));

    T exp_element =
        -half * (element_bb - ((element_ba * element_ba) / element_sum_aa_ob));

    return element_before_exp * exp(exp_element);
  }

  void reconstruction(int &iteration) {

    typename std::vector<event<T> >::iterator event_it;

    for (unsigned pixel_iterator = 0; pixel_iterator < n_pixels * n_pixels;
         ++pixel_iterator) {

      for (event_it = event_list.begin(); event_it != event_list.end();
           ++event_it) {

        T k = zero;
        T tan = event_tan(event_it->z_u, event_it->z_d);
        T y = event_y(event_it->dl, tan);
        T z = event_z(event_it->z_u, event_it->z_d, y, tan);
        /*
         *   loop -> event in bounding box(ellipse)
         *
         *
         *
         *
         */
      }
    }
  }

  T event_tan(T &z_u, T &z_d) const { return (z_u - z_d) / (2.0 * R_distance); }
  T event_y(T &dl, T &tan_event) const {
    return -half * (dl / sqrt(1 + (tan_event * tan_event)));
  }
  T event_z(T &z_u, T &z_d, T &y, T &tan_event) const {
    return half * (z_u + z_d + (two * y * tan_event));
  }

  pixel_location in_pixel(T &y, T &z) {

    return std::make_pair(((Scentilator_length + y) / pixel_size),
                          ((Scentilator_length + z) / pixel_size));
  }

  void ellipse_quadratic_form() {}
};

#endif  // SPET_RECONSTRUCTION_H
