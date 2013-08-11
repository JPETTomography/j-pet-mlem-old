#ifndef SPET_RECONSTRUCTION_H
#define SPET_RECONSTRUCTION_H

#include <cmath>
#include <vector>
#include <algorithm>
#include "event.h"

using std::fpclassify;

template <typename T = double> class Reconstruction {

  typedef std::pair<int, int> Pixel;
  typedef std::pair<T, T> Point;

private:
  static constexpr const T INVERSE_PI = T(1.0 / M_PI);
  static constexpr const T INVERSE_POW_TWO_PI = T(1.0 / (2.0 * M_PI * M_PI));

  T R_distance;
  T Scentilator_length;
  int n_pixels;
  T pixel_size;
  T sigma_z;
  T sigma_dl;
  T pow_sigma_z;
  T pow_sigma_dl;
  int m_pixel;
  T epsilon_distance;
  T half_pixel;

  std::vector<std::vector<event<T> > > event_list;
  std::vector<std::vector<T> > inverse_correlation_matrix;
  T sqrt_det_correlation_matrix;
  std::vector<T> rho;
  std::vector<T> rho_temp;
  std::vector<T> temp_kernels;

public:
  Reconstruction(T &R_distance, T &Scentilator_length, int &n_pixels,
                 T &pixel_size, T &sigma_z, T &sigma_dl)
      : R_distance(R_distance), Scentilator_length(Scentilator_length),
        n_pixels(n_pixels), pixel_size(pixel_size), sigma_z(sigma_z),
        sigma_dl(sigma_dl) {

    event_list.resize(n_pixels * n_pixels);
    rho.assign(n_pixels * n_pixels, T(1.0));
    rho_temp.assign(n_pixels * n_pixels, T(1.0));
    pow_sigma_z = sigma_z * sigma_z;
    pow_sigma_dl = sigma_dl * sigma_dl;
    m_pixel = (R_distance / pixel_size);
    epsilon_distance = pixel_size - T(0.000001);
    half_pixel = pixel_size / T(2.0);

#if DEBUG
    std::cout << "TEST ELIPSY: " << std::endl;

    T angle = T(32 * (M_PI / 180));

    point ellipse_center = std::pair<T, T>(T(0.0), T(0.0));
    point p = std::pair<T, T>(T(64.0), T(0));

    T _sin = std::sin(angle);
    T _cos = std::cos(angle);

    bool t = in_ellipse(_sin, _cos, ellipse_center, p);

    std::cout << "VALUE: " << t << std::endl;

    EllipseBoundingBox(ellipse_center, angle);

#endif

    inverse_correlation_matrix.resize(3, std::vector<T>(3, T()));

    inverse_correlation_matrix[0][0] = T(1) / pow_sigma_z;
    inverse_correlation_matrix[0][1] = T(0.0f);
    inverse_correlation_matrix[0][2] = T(0.0f);

    inverse_correlation_matrix[1][0] = T(0.0f);
    inverse_correlation_matrix[1][1] = T(1) / pow_sigma_z;
    inverse_correlation_matrix[1][2] = T(0.0f);

    inverse_correlation_matrix[2][0] = T(0.0f);
    inverse_correlation_matrix[2][1] = T(0.0f);
    inverse_correlation_matrix[2][2] = T(1) / pow_sigma_dl;

    sqrt_det_correlation_matrix = std::sqrt(inverse_correlation_matrix[0][0] *
                                            inverse_correlation_matrix[1][1] *
                                            inverse_correlation_matrix[2][2]);
  }

  T multiply_elements(std::vector<T> &vec_a, std::vector<T> &vec_b) {

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

  /*
  T multiply_elements(std::vector<T>& vec_a, std::vector<T>& vec_b) {

    std::vector<T> a(vec_a.size(), T(0));

    float output = 0.f;
    for (unsigned i = 0; i < 3; ++i) {
      a[i] += vec_a[i] * inverse_correlation_matrix[i][i];
      output += a[i] * vec_b[i];
    }

    return output;
  }
  */

  T kernel(T &y, T &z, T &angle, Point &pixel_center) {

    T _tan = std::tan(angle);
    T inv_cos = T(1) / std::cos(angle);
    T pow_inv_cos = inv_cos * inv_cos;

    std::vector<T> vec_o(3, T());
    std::vector<T> vec_a(3, T());
    std::vector<T> vec_b(3, T());

    vec_o[0] = -(pixel_center.first + y - R_distance) * _tan * pow_inv_cos;
    vec_o[1] = -(pixel_center.first + y + R_distance) * _tan * pow_inv_cos;
    vec_o[2] =
        -(pixel_center.first + y) * inv_cos * (T(1) + T(2) * (_tan * _tan));

    vec_a[0] = -(pixel_center.first + y - R_distance) * pow_inv_cos;
    vec_a[1] = -(pixel_center.first + y + R_distance) * pow_inv_cos;
    vec_a[2] = -T(2) * (pixel_center.first + y) * (inv_cos * _tan);

    vec_b[0] = pixel_center.second - (pixel_center.first * _tan);
    vec_b[1] = pixel_center.second - (pixel_center.first * _tan);
    vec_b[2] = -T(2) * pixel_center.first * inv_cos;

    T a_ic_a = multiply_elements(vec_a, vec_a);
    T b_ic_a = multiply_elements(vec_b, vec_a);
    T b_ic_b = multiply_elements(vec_b, vec_b);
    T o_ic_b = multiply_elements(vec_o, vec_b);

    T norm = a_ic_a + (T(2) * o_ic_b);

    T element_before_exp = INVERSE_POW_TWO_PI *
                           (sqrt_det_correlation_matrix / std::sqrt(norm)) *
                           sensitivity(y, z);

    T exp_element = -T(0.5) * (b_ic_b - ((b_ic_a * b_ic_a) / norm));

    T _exp = std::exp(exp_element);

    return (element_before_exp * _exp) /
           sensitivity(pixel_center.first, pixel_center.second);
  }

  T sensitivity(T y, T z) {

    T L_plus = (Scentilator_length / T(0.5) + z);
    T L_minus = (Scentilator_length / T(0.5) - z);
    T R_plus = R_distance + y;
    T R_minus = R_distance - y;

    return INVERSE_PI *
           (std::atan(std::max(-L_plus / R_minus, -L_plus / R_plus)) -
            std::atan(std::min(L_minus / R_minus, L_plus / R_plus)));
  }

  void operator()() {

    int k = 0;
    for (auto &col : event_list) {
      for (auto &row : col) {

        T tan = event_tan(row.z_u, row.z_d);
        T y = event_y(row.dl, tan);
        T z = event_z(row.z_u, row.z_d, y, tan);

        std::cout << "tan:= " << tan << " y:= " << y << " z:= " << z
                  << std::endl;

        Point ellipse_center = Point(y, z);
        T angle = std::atan(tan);
        bb_pixel_updates(ellipse_center, angle, y, z);
      }
    }
    rho = rho_temp;

    for (auto &it : rho) {

      if (it != T(1.0)) {

        std::cout << it << std::endl;
      }
    }
  }

  void test() {

    T y = 2.5;
    T z = 2.5;

    Point pixel = pixel_location(y, z);
    Point pp = pixel_center(pixel.first, pixel.second);

    std::cout << y << " " << z << std::endl;
    std::cout << "Pixel: " << pixel.first << " " << pixel.second << std::endl;
    std::cout << "P_center: " << pp.first << " " << pp.second << std::endl;

    Point ellipse_center = pixel_center(pixel.first, pixel.second);

    T angle = (0 * (M_PI / 180));

    // T _sin = std::sin(angle);
    // T _cos = std::cos(angle);

    T k = bb_pixel_updates(ellipse_center, angle, y, z);
    std::cout << "VALUE:" << k << std::endl;

    Pixel kk = pixel_location(200, 200);
    Point p = pixel_center(kk.first, kk.second);
    std::cout << "p: " << p.first << " " << p.second << std::endl;
    float ss = kernel(y, z, angle, p);

    if (fpclassify(ss) == FP_SUBNORMAL) {

      std::cout << "DUPA DENORMAL" << std::endl;
    }

    std::cout << "SS: " << ss << std::endl;
  }

  T bb_pixel_updates(Point &ellipse_center, T &angle, T y, T z) {

    std::cout << "el_cent in BB: " << ellipse_center.first << " "
              << ellipse_center.second << std::endl;

    T acc = T(0.0);

    T _cos = std::cos(angle);
    T _sin = std::sin(angle);

    T ux = sigma_z * _cos;
    T uy = sigma_z * _sin;
    T vx = sigma_dl * std::cos(angle + M_PI_2);
    T vy = sigma_dl * std::sin(angle + M_PI_2);

    T bbox_halfwidth = std::sqrt(ux * ux + vx * vx);
    T bbox_halfheight = std::sqrt(uy * uy + vy * vy);

    Pixel dl = pixel_location(ellipse_center.second - bbox_halfwidth,
                              ellipse_center.first - bbox_halfheight);
    Pixel ur = pixel_location(ellipse_center.second + bbox_halfwidth,
                              ellipse_center.first + bbox_halfheight);

    int iterator = 0;

    std::vector<std::pair<Pixel, T> > bb_pixels;

    for (int i = dl.first; i <= ur.first; ++i) {
      for (int j = dl.second; j <= ur.second; ++j) {

        Point p = pixel_center(j, i);

        int k = pixel_in_ellipse(p.first, p.second, ellipse_center, _sin, _cos);
        if (k == 2 || k == 1) {
          iterator++;

          T event_kernel = kernel(y, z, angle, p);

          bb_pixels.push_back(std::pair<Pixel, T>(Pixel(j, i), event_kernel));

          // std::cout << "event_kernel: " <<  event_kernel << " acc: " <<
          // event_kernel*sensitivity(p.first,p.second)*rho[mem_location(j,i)]
          // << std::endl;

          acc += event_kernel * sensitivity(p.first, p.second) *
                 rho[mem_location(j, i)];
        }
      }
    }

    for (auto &bb_event : bb_pixels) {

      rho_temp[mem_location(bb_event.first.first, bb_event.first.second)] =
          (bb_event.second *
           sensitivity(bb_event.first.first, bb_event.first.second)) /
          acc;

      //      std::cout << "BB_sum: " << (bb_event.second *
      // sensitivity(bb_event.first.first,bb_event.first.second))/acc <<
      // std::endl;
    }

    /*
    #if DEBUG_BBOX

        std::cout << "x_u: " << ellipse_center.second - bbox_halfwidth
                  << " y_u: " << ellipse_center.first - bbox_halfheight
                  << std::endl;
        std::cout << "x_d: " << ellipse_center.second + bbox_halfwidth
                  << " y_d: " << ellipse_center.first + bbox_halfheight
                  << std::endl;

    #endif
        Pixel dl = pixel_location(ellipse_center.second - bbox_halfwidth,
                                  ellipse_center.first - bbox_halfheight);
        Pixel ur = pixel_location(ellipse_center.second + bbox_halfwidth,
                                  ellipse_center.first + bbox_halfheight);

        //std::cout << "HERE!" << std::endl;
    #if DEBUG_BBOX
        std::cout << "dl : " << dl.first << " " << dl.second << std::endl;
        std::cout << "ur : " << ur.first << " " << ur.second << std::endl;
    #endif

    */

    return acc;
  }

  bool in_ellipse(T &_sin, T &_cos, Point &ellipse_center, Point &p) {

    T dy = (ellipse_center.first - p.first);
    T dz = (ellipse_center.second - p.second);
    T d1 = (_cos * dz + _sin * dy);
    T d2 = (_cos * dy + _sin * dz);

    // std::cout << "ELLIPSE VALUE: " << (d1 * d1 / pow_sigma_z) +
    //                                       (d2 * d2 / pow_sigma_dl) <<
    // std::endl;

    return ((d1 * d1 / pow_sigma_z) + (d2 * d2 / pow_sigma_dl)) <= T(1) ? true
                                                                        : false;
  }

  int pixel_in_ellipse(T &x, T &y, Point &ellipse_center, T &_sin, T &_cos) {

    Point pixel_ur = Point(x - half_pixel + epsilon_distance,
                           y - half_pixel + epsilon_distance);
    Point pixel_ul = Point(x - half_pixel, y - half_pixel + epsilon_distance);
    Point pixel_dr = Point(x - half_pixel + epsilon_distance, y - half_pixel);
    Point pixel_dl = Point(x - half_pixel, y - half_pixel);

#if DEBUG_PIXEL_IN_ELLIPSE
    std::cout << "Pixel_ur: " << pixel_ur.first << " " << pixel_ur.second
              << std::endl;
    std::cout << "Pixel_ul: " << pixel_ul.first << " " << pixel_ul.second
              << std::endl;
    std::cout << "Pixel_dr: " << pixel_dr.first << " " << pixel_dr.second
              << std::endl;
    std::cout << "Pixel_dl: " << pixel_dl.first << " " << pixel_dl.second
              << std::endl;
#endif

    bool ur = in_ellipse(_sin, _cos, ellipse_center, pixel_ur);
    bool ul = in_ellipse(_sin, _cos, ellipse_center, pixel_ul);
    bool dr = in_ellipse(_sin, _cos, ellipse_center, pixel_dr);
    bool dl = in_ellipse(_sin, _cos, ellipse_center, pixel_dl);

    // std::cout << ur << std::endl;
    // std::cout << ul << std::endl;
    // std::cout << dr << std::endl;
    // std::cout << dl << std::endl;

    if (ur && ul && dr && dl) {
      return 1;
    }
    if (ur || ul || dr || dl) {
      return 2;
    }

    return 0;
  }

  T event_tan(T &z_u, T &z_d) const {
    return (z_u - z_d) / (T(2) * R_distance);
  }
  T event_y(T &dl, T &tan_event) const {
    return -T(0.5) * (dl / sqrt(T(1) + (tan_event * tan_event)));
  }
  T event_z(T &z_u, T &z_d, T &y, T &tan_event) const {
    return T(0.5) * (z_u + z_d + (T(2) * y * tan_event));
  }

  // coord Plane
  Pixel pixel_location(T y, T z) {
    return Pixel(std::floor((R_distance + y) / pixel_size),
                 std::floor((R_distance + z) / pixel_size));
  }

  int mem_location(int y, int z) { return int(y * n_pixels + z); }

  // pixel Plane
  Point pixel_center(T y, T z) {
    return Point(
        (std::floor((y) * pixel_size - R_distance)) + (T(0.5) * pixel_size),
        (std::floor((z) * pixel_size - R_distance)) + (T(0.5) * pixel_size));
  }

  /* MIDPOINT ELLIPSE METHOD
  void DrawEllipse(int x0, int y0, int width, int height)
  {
      int a2 = width * width;
      int b2 = height * height;
      int fa2 = 4 * a2, fb2 = 4 * b2;
      int x, y, sigma;

      for (x = 0, y = height, sigma = 2*b2+a2*(1-2*height); b2*x <= a2*y; x++)
      {
        std::cout << x0 + x << " "  << y0 + y << std::endl;
        std::cout << x0 - x << " "  << y0 + y << std::endl;
        std::cout << x0 + x << " "  << y0 - y << std::endl;
        std::cout << x0 - x << " "  << y0 - y << std::endl;
          if (sigma >= 0)
          {
              sigma += fa2 * (1 - y);
              y--;
          }
          sigma += b2 * ((4 * x) + 6);
      }

      for (x = width, y = 0, sigma = 2*a2+b2*(1-2*width); a2*y <= b2*x; y++)
      {
        std::cout << x0 + x << " "  << y0 + y << std::endl;
        std::cout << x0 - x << " "  << y0 + y << std::endl;
        std::cout << x0 + x << " "  << y0 - y << std::endl;
        std::cout << x0 - x << " "  << y0 - y << std::endl;
          if (sigma >= 0)
          {
              sigma += fb2 * (1 - x);
              x--;
          }
          sigma += a2 * ((4 * y) + 6);
      }
  }
*/

  template <typename StreamType> Reconstruction &operator<<(StreamType &in) {

    event<T> temp_event;

    unsigned int n_pix;
    float pixel_s;
    unsigned int iter;
    unsigned int size;

    in >> n_pix;
    in >> pixel_s;
    in >> iter;
    in >> size;
    int i = 0;
    std::cout << "DATA:" << n_pix << " " << pixel_s << " " << iter << " "
              << size << std::endl;
#if DEBUG == 1
    std::cout << number_of_pixels << " " << pixel_s << " " << iter << " "
              << number_of_event_in_file << std::endl;
    std::cout << "VECTOR SIZE: " << event_list.size() << std::endl;

#endif
    for (;;) {

      if (in.eof()) {
        break;
      }
      i++;

      T z_u, z_d, dl;

      in >> z_u >> z_d >> dl;

      temp_event.z_u = z_u;
      temp_event.z_d = z_d;
      temp_event.dl = dl;

      T tan = event_tan(z_u, z_d);
      T y = event_y(z_d, tan);
      T z = event_z(z_u, z_d, y, tan);

      Pixel p = pixel_location(y, z);
#if DEBUG == 1

      std::cout << "pixel:" << p.first << " " << p.second << std::endl;
      event_list[p.first * n_pixels + p.second].push_back(temp_event);
      ++i;
      std::cout << "I:" << i << std::endl;
      std::cout << "DATA: " << z_u << "  " << z_d << " " << dl << std::endl;
      std::cout << "LOCATIONS: "
                << "y: " << y << " z: " << z << " tan: " << tan << std::endl;
      std::cout << "PIXEL: " << p.first << " " << p.second << " "
                << p.first *n_pixels + p.second << std::endl;
      std::cout << "N_PIXELS: " << number_of_pixels << std::endl;
// assert(p.first < number_of_pixels || p.second < number_of_pixels);
#endif
      event_list[mem_location(p.first, p.second)].push_back(temp_event);
    }

    return *this;
  }

  template <typename StreamType>
      friend StreamType &operator>>(StreamType &in, Reconstruction &r) {
    r << in;
    return in;
  }
};

#endif // SPET_RECONSTRUCTION_H
