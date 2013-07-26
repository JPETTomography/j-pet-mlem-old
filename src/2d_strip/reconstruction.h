#ifndef SPET_RECONSTRUCTION_
#define SPET_RECONSTRUCTION_H

#include <cmath>
#include <vector>
#include <algorithm>
#include "event.h"

template <typename T = double>
class spet_reconstruction {

  typedef std::pair<int, int> pixel_location;

 private:
  static constexpr const T INVERSE_PI = T(0.31830988);
  static constexpr const T INVERSE_TWO_PI = T(0.15915494);

  T R_distance;
  T Scentilator_length;
  int n_pixels;
  T pixel_size;
  T sigma_z;
  T sigma_dl;

  std::vector<std::vector<event<T>>> event_list;
  std::vector<std::vector<T>> inverse_correlation_matrix;
  T sqrt_det_correlation_matrix;
  std::vector<T> rho;

 public:
  spet_reconstruction(T& R_distance, T& Scentilator_length, int& n_pixels,
                      T& pixel_size, T& sigma_z, T& sigma_dl)
      : R_distance(R_distance),
        Scentilator_length(Scentilator_length),
        n_pixels(n_pixels),
        pixel_size(pixel_size),
        sigma_z(sigma_z),
        sigma_dl(sigma_dl)
        {

    event_list.resize(n_pixels * n_pixels);
    rho.resize(n_pixels * n_pixels);

    set_inverse_correlation_matrix();
  }

  void set_inverse_correlation_matrix() {

    inverse_correlation_matrix.resize(3, std::vector<T>(3, T()));

    T pow_sigma_z = sigma_z * sigma_z;
    T pow_sigma_dl = sigma_dl * sigma_dl;

    sqrt_det_correlation_matrix =
        std::sqrt(pow_sigma_z * pow_sigma_dl * pow_sigma_dl);

    inverse_correlation_matrix[0][0] =
        T(1) / pow_sigma_z;
    inverse_correlation_matrix[0][1] = T(0.0f);
    inverse_correlation_matrix[0][2] = T(0.0f);

    inverse_correlation_matrix[1][0] = T(0.0f);
    inverse_correlation_matrix[1][1] =
        T(1) / pow_sigma_z;
    inverse_correlation_matrix[1][2] = T(0.0f);

    inverse_correlation_matrix[2][0] = T(0.0f);
    inverse_correlation_matrix[2][1] = T(0.0f);
    inverse_correlation_matrix[2][2] =
        T(1) / pow_sigma_dl;
  }

  T multiply_elements(std::vector<T>& vec_a, std::vector<T>& vec_b) {

    std::vector<T> a(vec_a.size(), T(0));

    float output = 0.f;
    // add AVX
    for (unsigned i = 0; i < vec_a.size(); ++i) {
      for (unsigned j = 0; j < vec_b.size(); j++) {
        a[i] += vec_a[j] * inverse_correlation_matrix[j][i];
      }
      output += a[i] * vec_b[i];
    }

    return output;
  }

  T kernel(T& y, T& tan, pixel_location& pixel) {

    T angle = std::atan(tan);

    T event_inv_cos = T(1) / std::cos(angle * 180 * INVERSE_PI);
    T pow_event_inv_cos = event_inv_cos * event_inv_cos;

    std::vector<T> vec_o(3, T());
    std::vector<T> vec_a(3, T());
    std::vector<T> vec_b(3, T());

    vec_o[0] = -(pixel.first + y - R_distance) * tan * pow_event_inv_cos;
    vec_o[1] = -(pixel.first + y + R_distance) * tan * pow_event_inv_cos;
    vec_o[2] = -(pixel.first + y) * event_inv_cos * (T(1) + (T(2) * tan * tan));

    vec_a[0] = -(pixel.first + y - R_distance) * pow_event_inv_cos;
    vec_a[1] = -(pixel.first + y + R_distance) * pow_event_inv_cos;
    vec_a[2] = -(pixel.first + y) * event_inv_cos * tan;

    vec_b[0] = pixel.second - (pixel.first * tan);
    vec_b[1] = pixel.second - (pixel.first * tan);
    vec_b[2] = -T(2) * pixel.first * event_inv_cos;

    T element_aa = multiply_elements(vec_a, vec_a);
    T element_ba = multiply_elements(vec_b, vec_a);
    T element_bb = multiply_elements(vec_b, vec_b);
    T element_ob = multiply_elements(vec_o, vec_b);

    T element_sum_aa_ob = element_aa + (T(2) * element_ob);

    T element_before_exp = INVERSE_TWO_PI * (sqrt_det_correlation_matrix /
                                             std::sqrt(element_sum_aa_ob));

    T exp_element = -T(0.5) * (element_bb -
                               ((element_ba * element_ba) / element_sum_aa_ob));

    T _exp = std::exp(exp_element);

    return (element_before_exp *
            (sensitivity(pixel.first, pixel.second) * INVERSE_PI) * _exp);
  }

  T sensitivity(T y, T z) {

    pixel_location p = pixel_center(y, z);

    T L_plus = (Scentilator_length / T(0.5) + p.second);
    T L_minus = (Scentilator_length / T(0.5) - p.second);
    T R_plus = R_distance + p.first;
    T R_minus = R_distance - p.first;

    return INVERSE_PI *
           (std::atan(std::max(
                -L_plus / R_minus,
                ((-Scentilator_length) / T(0.5) + p.second) / R_plus)) -
            std::atan(std::min(L_minus / R_minus, L_plus / R_minus)));
  }

  void load_input(std::string fn) {

    ibstream in(fn, std::ios::binary);
    event<T> temp_event;

    unsigned int n_pix;
    float pixel_s;
    unsigned int iter;
    unsigned int size;

    in >> n_pix;
    in >> pixel_s;
    in >> iter;
    in >> size;

    std::cout <<"DATA:" <<  n_pix << " " << pixel_s << " " << iter << " "
              << size << std::endl;
#if DEBUG == 1
    std::cout << number_of_pixels << " " << pixel_s << " " << iter << " "
              << number_of_event_in_file << std::endl;
    std::cout << "VECTOR SIZE: " << event_list.size() << std::endl;

#endif
    for (;;) {

      T z_u, z_d, dl;

      in >> z_u >> z_d >> dl;
      std::cout << z_u << " " << z_d << " " << dl << std::endl;
      temp_event.z_u = z_u;
      temp_event.z_d = z_d;
      temp_event.dl = dl;

      T tan = get_event_tan(z_u, z_d);
      T y = get_event_y(z_d, tan);
      T z = get_event_z(z_u, z_d, y, tan);

      pixel_location p = in_pixel(y, z);

      std::cout << "pixel:" << p.first << " " << p.second << std::endl;

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

    std::cout << "VECTOR" << std::endl;

    typename std::vector<std::vector<event<T>>>::iterator it;
    typename std::vector<event<T>>::iterator jt;

    for (it = event_list.begin(); it != event_list.end(); ++it) {

      for (jt = (*it).begin(); jt != (*it).end(); ++jt) {

        std::cout << (*jt).z_u << " " << (*jt).z_d << " " << (*jt).dl
                  << std::endl;
      }
    }
  }

  void reconstruction(int& iteration) {
    std::cout << "START:" << std::endl;
    /*
    typename std::vector<std::vector<event<T>>>::iterator event_list_per_pixel;
    typename std::vector<event<T>>::iterator event_in_list;

    for (event_list_per_pixel = event_list; event_list_per_pixel != event_list.end();
         ++event_list_per_pixel) {

      for (event_in_list = *event_list_per_pixel.begin(); event_in_list != *event_list_per_pixel.end();
           ++event_in_list) {
    */
    for(auto& event_list_per_pixel : event_list){
      for(auto& event_per_pixel : event_list_per_pixel){

        std::cout << "HERE: " << std::endl;

        T k = T();
        T tan = get_event_tan(event_per_pixel.z_u, event_per_pixel.z_d);
        T y = get_event_y(event_per_pixel.dl, tan);
        T z = get_event_z(event_per_pixel.z_u, event_per_pixel.z_d, y, tan);
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

  T get_event_tan(T& z_u, T& z_d) const {
    return (z_u - z_d) / (T(2) * R_distance);
  }
  T get_event_y(T& dl, T& tan_event) const {
    return -T(0.5) * (dl / sqrt(T(1) + (tan_event * tan_event)));
  }
  T get_event_z(T& z_u, T& z_d, T& y, T& tan_event) const {
    return T(0.5) * (z_u + z_d + (T(2) * y * tan_event));
  }

  pixel_location in_pixel(T& y, T& z) {

    return std::make_pair(std::floor((R_distance - y) / pixel_size),
                          std::floor((R_distance - z) / pixel_size));
  }

  pixel_location pixel_center(T& y, T& z) {

    return std::make_pair(y * pixel_size + (T(0.5) * pixel_size),
                          z * pixel_size + (T(0.5) * pixel_size));
  }

  void ellipse_quadratic_form() {}
};

#endif  // SPET_RECONSTRUCTION_H
