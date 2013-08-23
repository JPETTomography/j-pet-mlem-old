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

  int iteration;
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

  std::vector<event<T>> event_list;
  std::vector<std::vector<T> > inverse_correlation_matrix;
  T sqrt_det_correlation_matrix;
  std::vector<std::vector<T>> rho;
  std::vector<std::vector<T>> rho_temp;
  std::vector<T> temp_kernels;

public:
  Reconstruction(int iteration, T R_distance, T Scentilator_length,
                 int n_pixels, T pixel_size, T sigma_z, T sigma_dl)
      : iteration(iteration), R_distance(R_distance),
        Scentilator_length(Scentilator_length), n_pixels(n_pixels),
        pixel_size(pixel_size), sigma_z(sigma_z), sigma_dl(sigma_dl) {

    rho.assign(n_pixels,std::vector<T>(n_pixels, T(10)));
    rho_temp.assign(n_pixels,std::vector<T>(n_pixels, T(10)));
    pow_sigma_z = sigma_z * sigma_z;
    pow_sigma_dl = sigma_dl * sigma_dl;
    m_pixel = (R_distance / pixel_size);
    epsilon_distance = pixel_size - T(0.000001);
    half_pixel = pixel_size / T(2.0);

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

  T kernel(T y, T z, T angle, Point pixel_center) {


      //std::cout << "INSIDE: " << y << " "  << z << " " << angle << " "  << pixel_center.first << " "  << pixel_center.second << std::endl;

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

    //std::cout << vec_o[0] << " " << vec_o[1] << " " << vec_o[2] << std::endl;
    vec_a[0] = -(pixel_center.first + y - R_distance) * pow_inv_cos;
    vec_a[1] = -(pixel_center.first + y + R_distance) * pow_inv_cos;
    vec_a[2] = -T(2) * (pixel_center.first + y) * (inv_cos * _tan);

    //std::cout << vec_a[0] << " " << vec_a[1] << " " << vec_a[2] << std::endl;

    vec_b[0] = pixel_center.second - (pixel_center.first * _tan);
    vec_b[1] = pixel_center.second - (pixel_center.first * _tan);
    vec_b[2] = -T(2) * pixel_center.first * inv_cos;

    //std::cout << vec_b[0] << " " << vec_b[1] << " " << vec_b[2] << std::endl;

    T a_ic_a = multiply_elements(vec_a, vec_a);
    T b_ic_a = multiply_elements(vec_b, vec_a);
    T b_ic_b = multiply_elements(vec_b, vec_b);
    T o_ic_b = multiply_elements(vec_o, vec_b);

    T norm = a_ic_a + (T(2) * o_ic_b);

    T element_before_exp = INVERSE_POW_TWO_PI *
                           (sqrt_det_correlation_matrix / std::sqrt(norm))* sensitivity(pixel_center.first,pixel_center.second);

    T exp_element = -T(0.5) * (b_ic_b - ((b_ic_a * b_ic_a) / norm));

    //std::cout << "BEFORE EXP: " << element_before_exp << std::endl;
    //std::cout << "EXP: " << exp_element << std::endl;

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

            (std::atan(std::min(L_minus / R_minus, L_plus / R_plus)) - std::atan(std::max(-L_plus / R_minus, -L_plus / R_plus)));
  }

  void operator()() {

      for (int i = 0; i < iteration; i++) {


                std::cout << "ITERATION: " <<i << std::endl;
                for (auto &event : event_list) {

                    T tan = event_tan(event.z_u, event.z_d);
                    T y = event_y(event.dl, tan);
                    T z = event_z(event.z_u, event.z_d, y, tan);

                    T angle = std::atan(tan);

                    Point ellipse_center = Point(y,z);

                    Pixel pp = pixel_location(y,z);

                    //rho_temp[pp.first][pp.second]+= T(1000);

                    bb_pixel_updates(ellipse_center, angle, y, z);
                }

      rho = rho_temp;

}
    std::ofstream file;
    file.open("pixels_output.txt");
    for (int x = 0; x < n_pixels; ++x) {
      for (int y = 0; y < n_pixels; ++y) {

        /*  if (rho[mem_location(y, x)] < 0.000001) {

            rho[mem_location(y, x)] = T(0.0);
          }
        */
          file << x << " " << y << " " << rho[x][y] << std::endl;
      }
    }

    // output reconstruction PNG
    png_writer png("out.png");
    png.write_header<>(n_pixels, n_pixels);

    T output_max = 0.0;
    for (auto &col : rho) {
      for (auto &row : col) {
        output_max = std::max(output_max, row);

     }
    }

    auto output_gain =
        static_cast<double>(std::numeric_limits<uint8_t>::max()) / output_max;

    for (int y = 0; y < n_pixels ; ++y) {
      uint8_t row[n_pixels];
      for (auto x = 0; x < n_pixels; ++x) {
        row[x] = std::numeric_limits<uint8_t>::max() -
                output_gain * rho[y][x];
      }
      png.write_row(row);
    }
  }

  void test() {
/*
      Point ellipse_center = Point(-19.8911,-109.594);

      Pixel p = pixel_location(ellipse_center.first,ellipse_center.second);

      std::cout << p.first << " " << p.second << std::endl;

   T angle = 0.192728;
   T _sin = std::sin(angle);
   T _cos = std::cos(angle);



   Point k = pixel_center(p.first,p.second);
   std::cout << "srodek: " << k.first << " " << k.second << std::endl;


   std::cout << "ELIPSA:" << in_ellipse(_sin, _cos, ellipse_center, Point(ellipse_center.first, ellipse_center.second)) << std::endl;

   bb_pixel_updates(ellipse_center,angle,ellipse_center.first,ellipse_center.second);

   std::cout << "PIXEL:" << p.first << " " << p.second << std::endl;
   std::cout << "KERNEL: " << kernel(ellipse_center.first,ellipse_center.second,angle,pixel_center(104,78)) << std::endl;
*/

 // bb_pixel_updates(ellipse_center,angle,ellipse_center.first,ellipse_center.second);

      std::cout << "---------------------------------------------" <<std::endl;

      Point pp = Point(0,0);

      Pixel p = pixel_location(pp.first,pp.second);

      Point pc = pixel_center(p.first,p.second);

      std::cout << "POINT: " << pp.first << " " << pp.second << " " << p.first << "  " << p.second << std::endl;

      std::cout << "CENTER: " << pc.first << " " << pc.second << std::endl;


      std::cout << "---------------------------------------------" <<std::endl;

      pp = Point(6.0,6.0);

      p = pixel_location(pp.first,pp.second);

      pc = pixel_center(p.first,p.second);

      std::cout << "POINT: " << pp.first << " " << pp.second << " " << p.first << "  " << p.second << std::endl;

      std::cout << "CENTER: " << pc.first << " " << pc.second << std::endl;

            std::cout << "---------------------------------------------" <<std::endl;

      pp = Point(-6.0,-6.0);

      p = pixel_location(pp.first,pp.second);

      pc = pixel_center(p.first,p.second);

      std::cout << "POINT: " << pp.first << " " << pp.second << " " << p.first << "  " << p.second << std::endl;

      std::cout << "CENTER: " << pc.first << " " << pc.second << std::endl;

            std::cout << "---------------------------------------------" <<std::endl;


      pp = Point(-6.0,6.0);

      p = pixel_location(pp.first,pp.second);

      pc = pixel_center(p.first,p.second);

      std::cout << "POINT: " << pp.first << " " << pp.second << " " << p.first << "  " << p.second << std::endl;

      std::cout << "CENTER: " << pc.first << " " << pc.second << std::endl;

            std::cout << "---------------------------------------------" <<std::endl;

      pp = Point(6.0,-6.0);

      p = pixel_location(pp.first,pp.second);

      pc = pixel_center(p.first,p.second);

      std::cout << "POINT: " << pp.first << " " << pp.second << " " << p.first << "  " << p.second << std::endl;

      std::cout << "CENTER: " << pc.first << " " << pc.second << std::endl;
  }

  void bb_pixel_updates(Point &ellipse_center, T &angle, T y, T z) {

      T acc = T(0.0);

      T _sin = std::sin(angle);
      T _cos = std::cos(angle);

      Pixel center_pixel = pixel_location(ellipse_center.first,ellipse_center.second);

      //std::cout << "CENTER PIXEL: " << center_pixel.first << " " << center_pixel.second << std::endl;

      Pixel ur = Pixel(center_pixel.first - 25,center_pixel.second + 25);
      Pixel dl = Pixel(center_pixel.first + 25,center_pixel.second - 25);


      //std::cout << ur.first << " " << ur.second << " " << dl.first << " " << dl.second << std::endl;

    //std::cout << "LIMIT: " << dl.second << " " << ur.second << std::endl;
     //std::cout << "START:" << ur.first << " " << dl.first << std::endl;


      std::vector<std::pair<Pixel,T>> ellipse_kernels;

      for(int iz = dl.second; iz < ur.second;++iz ){
          for(int iy = ur.first; iy < dl.first;++iy){

              Point  pp = pixel_center(iy,iz);

              if(in_ellipse(_sin,_cos,ellipse_center,pp)){


               //std::cout<<"PIXEL: " << iy << " " << iz << std::endl;


                  T event_kernel = kernel(y,z,angle,pp);


                  ellipse_kernels.push_back(std::pair<Pixel,T>(Pixel(iy,iz),event_kernel));
                 // std::cout << "EVENT: " << y << " " << z << " ANGLE: " << angle <<  " EVENT PIXEL: " << center_pixel.first << " " << center_pixel.second << std::endl;
                 // std::cout << "PIXEL: " << iy << " " << iz << " PIXEL CENTER: " << pp.first << "  "  << pp.second << std::endl;
                 // std::cout << "Kernel:  " << event_kernel << " sensitivity:  " << sensitivity(pp.first,pp.second) << std::endl;

                  acc+= event_kernel*sensitivity(pp.first,pp.second)*rho[iy][iz];

                 //std::cout << "ACC: " << acc << std::endl;



              }

          }

      }

      for(auto& e: ellipse_kernels){

        if(e.second > 0){rho_temp[e.first.first][e.first.second] += e.second*rho[e.first.first][e.first.second]/acc;
}
         // std::cout << e.second*rho[e.first.first][e.first.second]/acc << std::endl;

      }

 }

  bool in_ellipse(T _sin, T _cos, Point ellipse_center, Point p) {

    T dy = (p.first - ellipse_center.first);
    T dz = (p.second - ellipse_center.second);
    T d1 = (_cos * dz + _sin * dy);
    T d2 = (_sin * dz - _cos * dy);

    T tg = (_sin / _cos);

    T A = (((T(4.0) * (T(1.0) / (_cos * _cos))) / pow_sigma_dl) +
            (T(2.0) * tg * tg / pow_sigma_z));
    T B = -T(4.0) * tg / pow_sigma_z;
    T C = T(2.0) / pow_sigma_z;

      return (((A * (d2 * d2)) + (B * d1 * d2) + (C * (d1 * d1)))) <= T(1) ? true
                                                                      : false;
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
    return Pixel(std::floor((R_distance - y) / pixel_size),
                 std::ceil((R_distance + z) / pixel_size));
  }

  int mem_location(int y, int z) { return int(y * n_pixels + z); }

  // pixel Plane
  Point pixel_center(T y, T z) {

    int sgn_y =  sgn<T>(y);
    int sgn_z =  sgn<T>(z);


    return Point(
        (std::floor(R_distance - (y) * pixel_size)) - (sgn<T>(y)*T(0.5) * pixel_size),
        (std::floor((z) * pixel_size - R_distance)) - (sgn<T>(z)*T(0.5) * pixel_size));
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
      event_list.push_back(temp_event);
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
