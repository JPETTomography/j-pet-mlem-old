#ifndef SPET_RECONSTRUCTION_
#define SPET_RECONSTRUCTION_H

#include <cmath>
#include <vector>
#include <algorithm>
#include "data_structures.h"

template <typename T = float> class spet_reconstruction {

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
  std::vector<std::vector<event<T>>> event_list;
  std::vector<std::vector<T>> inverse_correlation_matrix;
  T inverse_correlation_matrix_factor;
  std::vector<T> rho;
  static constexpr const T INVERSE_PI = T(0.31830988);
  static constexpr const T INVERSE_TWO_PI = T(0.15915494);

 public:
  spet_reconstruction(T& R_distance,
                      T& Scentilator_length,
                      int& n_pixels,
                      T& pixel_size,
                      T& sigma_z,
                      T& sigma_dl,
                      T& gamma)
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
  }

  void set_inverse_correlation_matrix() {

    inverse_correlation_matrix.resize(3, std::vector<T>(3, T()));

    T pow_sigma_z = sigma_z * sigma_z;
    T pow_sigma_dl = sigma_dl * sigma_dl;
    T pow_gamma = gamma * gamma;

    std::cout << "pow_gammas: " << pow_sigma_z << " " << pow_sigma_dl
              << std::endl;

    inverse_correlation_matrix[0][0] =
        T(1) / pow_sigma_z;  // pow_gamma - (pow_sigma_z * pow_sigma_dl);
    inverse_correlation_matrix[0][1] = T();
    inverse_correlation_matrix[0][2] = T();

    std::cout << "D1: " << 1 / pow_sigma_z << std::endl;

    inverse_correlation_matrix[1][0] = T();
    inverse_correlation_matrix[1][1] =
        T(1) / pow_sigma_z;  // pow_gamma - (pow_sigma_z * pow_sigma_dl);
    inverse_correlation_matrix[1][2] = T();

    std::cout << "D2: " << T(1) / pow_sigma_z << std::endl;

    inverse_correlation_matrix[2][0] = T();
    inverse_correlation_matrix[2][1] = T();
    inverse_correlation_matrix[2][2] =
        T(1) / pow_sigma_dl;  // -pow_sigma_z * pow_sigma_z;

    std::cout << "D3: " << T(1) / pow_sigma_dl << std::endl;

    inverse_correlation_matrix_factor =
        T(1) / ((pow_gamma * pow_sigma_z) + (pow_gamma * pow_sigma_z) +
                (pow_sigma_z * pow_sigma_z * pow_sigma_dl));
  }

  T multiply_elements(std::vector<T>& vec_a, std::vector<T>& vec_b) {

    T a[vec_a.size()];
    T output = T();
    // add AVX
    for (unsigned i = 0; i < vec_a.size(); ++i) {
      for (unsigned j = 0; j < vec_a.size(); j++) {

        a[i] += vec_a[j] * inverse_correlation_matrix[j][i];
        std::cout << "a[i]: " << a[i] << std::endl;
      }
      output += a[i] * vec_b[i];
    }

    return output;
  }

  T kernel(T& y, T& tan, pixel_location& pixel) {

    std::cout << pixel.first << " " << pixel.second << std::endl;

    T angle = std::atan(tan);

    T event_inv_cos = T(1) / std::cos(angle * 180 * INVERSE_PI);
    T pow_event_inv_cos = event_inv_cos * event_inv_cos;

    std::cout << "INVERSE_CORRELATION_MATIRX: "
              << inverse_correlation_matrix[0][0] << " "
              << inverse_correlation_matrix[1][1] << " "
              << inverse_correlation_matrix[2][2] << std::endl;

    std::cout << "INVERSE_CORRELATION_MATIX_FACTOR: "
              << inverse_correlation_matrix_factor << std::endl;

    std::cout << "angle: " << std::atan(tan) << " " << angle * 180 * INVERSE_PI
              << " inv_cos: " << event_inv_cos
              << " pow_inv_cos: " << pow_event_inv_cos << std::endl;

    std::vector<T> vec_o(3, T());
    std::vector<T> vec_a(3, T());
    std::vector<T> vec_b(3, T());

    vec_o[0] = -(pixel.first + y - R_distance) * tan * pow_event_inv_cos;
    vec_o[1] = -(pixel.first + y + R_distance) * tan * pow_event_inv_cos;
    vec_o[2] = -(pixel.first + y) * event_inv_cos * (T(1) + (T(2) * tan * tan));

    std::cout << "VEC_o: " << vec_o[0] << " " << vec_o[1] << " " << vec_o[2]
              << std::endl;

    vec_a[0] = -(pixel.first + y - R_distance) * pow_event_inv_cos;
    vec_a[1] = -(pixel.first + y + R_distance) * pow_event_inv_cos;
    vec_a[2] = -(pixel.first + y) * event_inv_cos * tan;

    std::cout << "VEC_a: " << vec_a[0] << " " << vec_a[1] << " " << vec_a[2]
              << std::endl;

    vec_b[0] = pixel.second - (pixel.first * tan);
    vec_b[1] = pixel.second - (pixel.first * tan);
    vec_b[2] = -T(2) * pixel.first * event_inv_cos;

    std::cout << "VEC_b: " << vec_b[0] << " " << vec_b[1] << " " << vec_b[2]
              << std::endl;

    T element_aa = multiply_elements(vec_a, vec_a);
    T element_ba = multiply_elements(vec_b, vec_a);
    T element_bb = multiply_elements(vec_b, vec_b);
    T element_ob = multiply_elements(vec_o, vec_b);

    std::cout << "element_aa: " << element_aa << " element_ob: " << element_ob
              << std::endl;

    T element_sum_aa_ob = element_aa + (T(2) * element_ob);

    T det_inverse_matrix = inverse_correlation_matrix[0][0] *
                           inverse_correlation_matrix[1][1] *
                           inverse_correlation_matrix[2][2];

    std::cout << "DET MATRIX: " << det_inverse_matrix << std::endl;

    T element_before_exp =
        INVERSE_TWO_PI * (det_inverse_matrix / sqrt(element_sum_aa_ob));

    std::cout << "ELEMENT_BEFORE_EXP: " << element_before_exp << std::endl;

    T exp_element =
        -T(0.5) * (element_bb - ((element_ba * element_ba) / element_sum_aa_ob));

    return element_before_exp * exp_element;
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
      std::cout << "LOCATIONS: "
                << "y: " << y << " z: " << z << " tan: " << tan << std::endl;
      std::cout << "PIXEL: " << p.first << " " << p.second << " "
                << p.first* n_pixels + p.second << std::endl;
      std::cout << "N_PIXELS: " << number_of_pixels << std::endl;
// assert(p.first < number_of_pixels || p.second < number_of_pixels);
#endif
      // event_list[p.first * n_pixels + p.second].push_back(temp_event);

      if (in.eof()) {
        break;
      }
    }

#if DEBUG == 1

    std::cout << "VECTOR" << std::endl;

    typename std::vector<std::vector<event<T>>>::iterator it;
    typename std::vector<event<T>>::iterator jt;

    for (it = event_list.begin(); it != event_list.end(); ++it) {

      for (jt = (*it).begin(); jt != (*it).end(); ++jt) {

        std::cout << (*jt).z_u << " " << (*jt).z_d << " " << (*jt).dl
                  << std::endl;
      }
    }
#endif
  }

  T sensitivity(T& y, T& z) {

    T a = (T(0.5) * Scentilator_length - z) / (R_distance - y);
    T b = (T(0.5) * Scentilator_length + z) / (R_distance + y);
    T c = (T(0.5) * Scentilator_length + z) / (R_distance - y);

    return INVERSE_PI *
           (std::atan(std::min(a, b)) - std::atan(std::min(-c, -b)));
  }

  void reconstruction(int& iteration) {

    typename std::vector<event<T>>::iterator event_it;

    for (unsigned pixel_iterator = 0; pixel_iterator < n_pixels * n_pixels;
         ++pixel_iterator) {

      for (event_it = event_list.begin(); event_it != event_list.end();
           ++event_it) {

        T k = T();
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

  T event_tan(T& z_u, T& z_d) const {
    return (z_u - z_d) / (T(2) * R_distance);
  }
  T event_y(T& dl, T& tan_event) const {
    return -T(0.5) * (dl / sqrt(T(1) + (tan_event * tan_event)));
  }
  T event_z(T& z_u, T& z_d, T& y, T& tan_event) const {
    return T(0.5) * (z_u + z_d + (T(2) * y * tan_event));
  }

  pixel_location in_pixel(T& y, T& z) {

    return std::make_pair((((Scentilator_length / T(2)) + y) / pixel_size),
                          (((Scentilator_length / T(2)) + z) / pixel_size));
  }

  void ellipse_quadratic_form() {}
};

#endif  // SPET_RECONSTRUCTION_H
