#ifndef SPET_RECONSTRUCTION_
#define SPET_RECONSTRUCTION_H

#include <cmath>
#include <vector>
#include <algorithm>
#include "event.h"

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

public:
  Reconstruction(T &R_distance, T &Scentilator_length, int &n_pixels,
                 T &pixel_size, T &sigma_z, T &sigma_dl)
      : R_distance(R_distance), Scentilator_length(Scentilator_length),
        n_pixels(n_pixels), pixel_size(pixel_size), sigma_z(sigma_z),
        sigma_dl(sigma_dl) {

    event_list.resize(n_pixels * n_pixels);
    rho.assign(n_pixels * n_pixels, T(0.1));
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

#if DEBUG_KERNEL
    std::cout << "ANGLE : " << angle << std::endl;
#endif

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

#if DEBUG_KERNEL
    std::cout << "KERNEL METHOD" << std::endl;
    std::cout << "y: " << y << " z: " << z
              << " p_center: " << pixel_center.first << " "
              << pixel_center.second << std::endl;
    std::cout << "ANGLE : " << angle << std::endl;
    std::cout << "vec_o: " << vec_o[0] << " " << vec_o[1] << " " << vec_o[2]
              << std::endl;
#endif

    vec_a[0] = -(pixel_center.first + y - R_distance) * pow_inv_cos;
    vec_a[1] = -(pixel_center.first + y + R_distance) * pow_inv_cos;
    vec_a[2] = -T(2) * (pixel_center.first + y) * (inv_cos * _tan);

#if DEBUG_KERNEL
    std::cout << "vec_a: " << vec_a[0] << " " << vec_a[1] << " " << vec_a[2]
              << std::endl;
#endif

    vec_b[0] = pixel_center.second - (pixel_center.first * _tan);
    vec_b[1] = pixel_center.second - (pixel_center.first * _tan);
    vec_b[2] = -T(2) * pixel_center.first * inv_cos;

#if DEBUG_KERNEL
    std::cout << "vec_b: " << vec_b[0] << " " << vec_b[1] << " " << vec_b[2]
              << std::endl;
#endif

    T a_ic_a = multiply_elements(vec_a, vec_a);
    T b_ic_a = multiply_elements(vec_b, vec_a);
    T b_ic_b = multiply_elements(vec_b, vec_b);
    T o_ic_b = multiply_elements(vec_o, vec_b);

    T norm = a_ic_a + (T(2) * o_ic_b);

    T element_before_exp =
        INVERSE_POW_TWO_PI * (sqrt_det_correlation_matrix / std::sqrt(norm)) *
        sensitivity(y, z);

#if DEBUG_KERNEL
    std::cout << "SQRT DET_COR := " << sqrt_det_correlation_matrix << std::endl;
    std::cout << "sqrt norm := " << std::sqrt(norm) << std::endl;
    std::cout << "1/2pi^2 * (sqrt((det)/norm)) " << element_before_exp
              << std::endl;
    std::cout << "ELEMENT BEFORE EXP: " << element_before_exp << std::endl;
#endif

    T exp_element = -T(0.5) * (b_ic_b - ((b_ic_a * b_ic_a) / norm));

    T _exp = std::exp(exp_element);

#if DEBUG_KERNEL
    std::cout << "b.ic.b := " << b_ic_b << std::endl;
    std::cout << "b.ic.a := " << b_ic_a << std::endl;
    std::cout << "a.ic.a := " << a_ic_a << std::endl;
    std::cout << "o.ic.b := " << o_ic_b << std::endl;
    std::cout << "norm: := " << norm << std::endl;
    std::cout << "0.5*((b.ic.b - (b.ic.a)^2/(a.ic.a + 2 o.ic.b ))):"
              << 0.5 * (b_ic_b - ((b_ic_a * b_ic_a) / norm)) << std::endl;

    std::cout << "x,y SENS: " << sensitivity(y, z) << std::endl;

    std::cout << "PIXEL SENS: " << sensitivity(pixel_center.first,
                                               pixel_center.second)
              << std::endl;

    std::cout << "EXP: " << _exp << std::endl;
#endif

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
        std::cout << "k:= " << k << std::endl;
        T tan = event_tan(row.z_u, row.z_d);
        T y = event_y(row.dl, tan);
        T z = event_z(row.z_u, row.z_d, y, tan);

        Point pcenter = pixel_center(y, z);

        //T main_kernel = kernel(y, z, tan, pcenter);
        //std::cout << "MAIN_KERNEL: " << main_kernel << std::endl;

        Point ellipse_center = std::pair<T, T>(y, z);
        T angle = std::atan(tan);
        T denominator = ellipse_bounding_box(ellipse_center, angle);
        ++k;
      }
    }
  }

  void test() {

    T x = 0.0;
    T y = 0.0;
    /*
    typename std::vector<std::vector<event<T>>>::iterator it;
    typename std::vector<event<T>>::iterator jt;
    int k = 0;
    for(it = event_list.begin(); it != event_list.end();++it){
        //for(jt = it->begin();jt != it->end();++jt){


            k+= (*it).size();

       // }

        std::cout << "k: " << k  << std::endl;



    }

    */

    Point pixel = pixel_location(x, y);
    Point pp = pixel_center(x, y);

    std::cout << "Pixel: " << pixel.first << " " << pixel.second << std::endl;
    std::cout << "P_center: " << pp.first << " " << pp.second << std::endl;

    Point ellipse_center = pixel_center(x, y);
    T angle = (0 * (M_PI / 180));

    //T _sin = std::sin(angle);
    //T _cos = std::cos(angle);

    T k = ellipse_bounding_box(ellipse_center, angle);
    std::cout << "VALUE:" << k << std::endl;

    //T a = 0;
    // T b = 0;

    //T ss = kernel(a,b,angle,pp);

  }

  T ellipse_bounding_box(Point &ellipse_center, T &angle) {

    T acc = T(0.0);

    T _cos = std::cos(angle);
    T _sin = std::sin(angle);

    T ux = sigma_z * _cos;
    T uy = sigma_z * _sin;
    T vx = sigma_dl * std::cos(angle + M_PI_2);
    T vy = sigma_dl * std::sin(angle + M_PI_2);

    T bbox_halfwidth = std::sqrt(ux * ux + vx * vx);
    T bbox_halfheight = std::sqrt(uy * uy + vy * vy);

    std::cout << "x_u: " << ellipse_center.second - bbox_halfwidth
              << " y_u: " << ellipse_center.first - bbox_halfheight
              << std::endl;
    std::cout << "x_d: " << ellipse_center.second + bbox_halfwidth
              << " y_d: " << ellipse_center.first + bbox_halfheight
              << std::endl;

    Pixel dl = pixel_location(ellipse_center.second - bbox_halfwidth,
                              ellipse_center.first - bbox_halfheight);
    Pixel ur = pixel_location(ellipse_center.second + bbox_halfwidth,
                              ellipse_center.first + bbox_halfheight);

    //std::cout << "HERE!" << std::endl;

    std::cout << "dl : " << dl.first << " " << dl.second << std::endl;
    std::cout << "ur : " << ur.first << " " << ur.second << std::endl;

    int iterator = 0;
    for (int i = dl.first; i <= ur.first; ++i) {
      for (int j = dl.second; j <= ur.second; ++j) {

        //   std::cout << "ELLIPSE CENTER: " << ellipse_center.first << " " <<
        // ellipse_center.second << std::endl;
        //   std::cout << "dl : " << dl.first << " " << dl.second << std::endl;
        //   std::cout << "ur : " << ur.first << " " << ur.second << std::endl;

        int location = mem_location(i, j);

        T px = (i * pixel_size) - (R_distance);
        T py = (j * pixel_size) - (R_distance);

        //  std::cout << "PIXEL (i,j): " << i << " " << j << std::endl;
        //std::cout << "PIXEL (x,y): " << px << " " << py << std::endl;

        Point pixel = pixel_center(px, py);

        //  std::cout << "PIXEL center (x,y): " << pixel.first << " " <<
        // pixel.second << std::endl;
        //  std::cout << "R_distance: " << R_distance << std::endl;
        //  std::cout << "Scentilator_length: " << Scentilator_length <<
        // std::endl;

        int k = pixel_in_ellipse(pixel.first, pixel.second, ellipse_center,
                                 _sin, _cos);

        // if (k != 0) {
        //  std::cout << "k: " << k << std::endl;
        // }

        if (k == 2 || k == 1) {
          std::cout << "P: " << i << " " << j << std::endl;
          iterator++;

          if (event_list[location].size() != 0) {
            for (unsigned int in = 0; in < event_list[location].size(); ++in) {

              T _tan = event_tan(event_list[location][in].z_u,
                                 event_list[location][in].z_d);
              T y = event_y(event_list[location][in].dl, _tan);
              T z = event_z(event_list[location][in].z_u,
                            event_list[location][in].z_d, y, _tan);

              std::cout << "INSIDE: " << _tan << " " << y << " " << z
                        << std::endl;

              Point p = pixel_center(y, z);
              T angle = std::atan(_tan);
              T kk = kernel(y, z, angle, p);
              std::cout << "kernel: " << kernel(y, z, angle, p) << std::endl;
              acc += kk;

            }
            //std::cout << iterator << std::endl;
          }

        }
      }
    }
    std::cout << "CALOSC: " << acc << std::endl;
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

    // std::cout << "Pixel_ur: " << pixel_ur.first << " " << pixel_ur.second <<
    // std::endl;
    // std::cout << "Pixel_ul: " << pixel_ul.first << " " << pixel_ul.second <<
    // std::endl;
    // std::cout << "Pixel_dr: " << pixel_dr.first << " " << pixel_dr.second <<
    // std::endl;
    // std::cout << "Pixel_dl: " << pixel_dl.first << " " << pixel_dl.second <<
    // std::endl;

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

  Pixel pixel_location(T y, T z) {
    return Pixel(std::floor((R_distance + y) / pixel_size),
                 std::floor((R_distance + z) / pixel_size));
  }

  int mem_location(int &y, int &z) { return int(y * n_pixels + z); }

  Point pixel_center(T y, T z) {
    return Point((std::floor((y) / pixel_size)) + (T(0.5) * pixel_size),
                 (std::floor((z) / pixel_size)) + (T(0.5) * pixel_size));
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
