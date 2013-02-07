#pragma once

#include <iostream>
#include <vector>
#include <list>

#include <ctime>

#include "detector_ring.h"
#include "sparse_matrix.h"
#include "pixel.h"

template <typename FType = double, typename SType = int> class Reconstruction {
 public:
  typedef FType F;
  typedef SType S;
  typedef std::pair<S, S> PixelLocation;
  typedef std::pair<S, S> LOR;
  typedef std::vector<F> Output;

  struct HitsPerPixel {
    S index;
    F probability;
  };

  struct MeanPerLOR {
    LOR lor;
    S n;
  };

  typedef uint32_t FileInt;
  typedef uint16_t FileHalf;

  Reconstruction(S n_iterations,
                 std::string matrix,
                 std::string mean,
                 F threshold = (F) 0.0)
      : n_iterations_(n_iterations),
        threshold_(threshold) {
    ibstream in(matrix);

    if (!in) {
      std::cerr << "error opening matrix file '" << matrix << "'" << std::endl;
      exit(-1);
    }

    typedef uint32_t FileInt;
    typedef uint16_t FileHalf;

    FileInt in_magic;

    in >> in_magic;

    if (in_magic != SparseMatrix<Pixel<>, LOR>::MAGIC_VERSION_FULL) {
      throw("invalid input system matrix file");
    }

    in >> n_pixels_;
    in >> emissions_;
    in >> n_tubes_;

    total_n_pixels_ = n_pixels_ * n_pixels_;
    rho_.resize(total_n_pixels_, (F) 0.0);
    rho_detected_.resize(total_n_pixels_, (F) 0.0);

    scale.resize(total_n_pixels_, (F) 0.0);

    std::ifstream mean_file(mean);
    if (!mean_file) {
      std::cerr << "error opening mean file '" << mean << "'" << std::endl;
      exit(-1);
    }
    std::vector<MeanPerLOR> lor_mean;

    clock_t start = clock();
    for (;;) {

      S x, y, value;

      mean_file >> y >> x >> value;
      if (mean_file.eof())
        break;

      MeanPerLOR temp_obj;

      LOR lor(x, y);

      temp_obj.lor = lor;
      temp_obj.n = value;

      lor_mean.push_back(temp_obj);
    }

    clock_t stop = clock();

    double time = static_cast<double>(stop - start) / CLOCKS_PER_SEC;
    std::cout << "means read time = " << time << "s" << std::endl;

    std::vector<HitsPerPixel> pixels;

    S index = 0;
    n_non_zero_elements_ = 0;

    start = clock();
    for (;;) {

      FileHalf a, b;
      in >> a >> b;
      if (in.eof())
        break;

      LOR lor(a, b);

      n.push_back(get_mean_per_lor(a, b, lor_mean));

      FileInt count;

      in >> count;

      for (S i = 0; i < count; ++i) {

        FileHalf x, y;
        FileInt hits;

        in >> x >> y >> hits;
        HitsPerPixel data;
        data.probability = static_cast<F>(hits / static_cast<F>(emissions_));
        data.index = location(x, y, n_pixels_);
        if (threshold_ > (F) 0.0) {
          if (data.probability < threshold_)
            data.probability = (F) 0.0;
          else
            data.probability = (F) 1.0;
        }
        scale[data.index] += data.probability;
        n_non_zero_elements_++;
        pixels.push_back(data);
      }

      system_matrix.push_back(pixels);
      pixels.clear();
      index++;
    }
    stop = clock();

    time = static_cast<double>(stop - start) / CLOCKS_PER_SEC;
    std::cout << "matrix read time = " << time << "s\n";

    for (auto it_vector = system_matrix.begin();
         it_vector != system_matrix.end(); it_vector++) {
      for (auto it_list = it_vector->begin(); it_list != it_vector->end();
           it_list++) {
        S pixel = it_list->index;
        it_list->probability /= scale[pixel];
      }
    }

    for (S p = 0; p < n_pixels_ * n_pixels_; ++p) {
      if (scale[p] > 0)
        rho_detected_[p] = (F) 1.0;
    }
    std::cout << "   Pixels: " << n_pixels_ << std::endl;
    std::cout << "Emissions: " << emissions_ << std::endl;
    std::cout << "Detectors: " << n_tubes_ << std::endl;
    std::cout << "     LORs: " << system_matrix.size() << std::endl;
    std::cout << "Non zero elements: " << n_non_zero_elements_ << std::endl;
  }

  ~Reconstruction() {}

  void emt(S n_iterations) {

    F y[n_pixels_* n_pixels_];
    F u;

    clock_t start = clock();

    for (S i = 0; i < n_iterations; ++i) {
      std::cout << ".";
      std::cout.flush();

      for (S p = 0; p < n_pixels_ * n_pixels_; ++p) {
        y[p] = (F) 0.0;
      }

      S t = 0;
      for (auto it_vector = system_matrix.begin();
           it_vector != system_matrix.end(); it_vector++) {
        u = (F) 0.0;
        if (n[t] > 0) {
          for (auto it_list = it_vector->begin(); it_list != it_vector->end();
               it_list++) {
            u += rho_detected_[it_list->index] * it_list->probability;
          }
          F phi = n[t] / u;
          for (auto it_list = it_vector->begin(); it_list != it_vector->end();
               ++it_list) {
            y[it_list->index] += phi * it_list->probability;
          }
        }
        t++;
      }

      for (S p = 0; p < n_pixels_ * n_pixels_; ++p) {
        if (scale[p] > 0) {
          rho_detected_[p] *= y[p];
        }
      }
    }
    clock_t stop = clock();
    std::cout << std::endl;

    for (S p = 0; p < n_pixels_ * n_pixels_; ++p) {
      if (scale[p] > 0) {
        rho_[p] = (rho_detected_[p] / scale[p]);
      }
    }

    double time = static_cast<double>(stop - start) / CLOCKS_PER_SEC;
    std::cout << "time = " << time << "s "
              << "time/iter = " << time / n_iterations << "s" << std::endl;
    std::cout << "op/sec = " << n_non_zero_elements_ * (n_iterations / time)
              << std::endl;
  }

  S get_n_pixels() { return n_pixels_; }

  F rho(S p) const { return rho_[p]; }
  F rho_detected(S p) const { return rho_detected_[p]; }
  std::vector<F> rho() const { return rho_; }
  std::vector<F> rho_detected() const { return rho_detected_; }

  void set_threshold(F t) { threshold_ = t; }

  static S location(S x, S y, S size) { return y * size + x; }

 private:
  S get_mean_per_lor(FileHalf& a, FileHalf& b, std::vector<MeanPerLOR>& mean) {
    for (auto it = mean.begin(); it != mean.end(); ++it) {
      if (a == it->lor.first && b == it->lor.second) {
        return it->n;
      }
    }

    return 0;
  }

  S n_tubes_;
  S n_pixels_;
  S total_n_pixels_;
  S n_iterations_;
  S emissions_;
  S n_non_zero_elements_;
  std::vector<std::vector<HitsPerPixel>> system_matrix;
  std::vector<S> list_of_lors;
  std::vector<S> n;
  std::vector<F> scale;

  std::vector<F> rho_;
  std::vector<F> rho_detected_;

  F threshold_;

};
