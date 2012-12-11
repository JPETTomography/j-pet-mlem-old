#pragma once


#include <iostream>
#include <vector>
#include <list>

#include<ctime>

#include "detector_ring.h"

inline int  LOCATION(int x,int y,int size)  {return y*size + x;}

template <typename F = double>
class reconstruction {
public:
  typedef std::pair<size_t, size_t> pixel_location;
  typedef std::pair<size_t, size_t> lor_type;
  typedef std::vector<F> output_type;

  struct hits_per_pixel {
    int index;
    F probability;
  };

  struct mean_per_lor {
    lor_type lor;
    int n;
  };

  typedef uint32_t file_int;
  typedef uint16_t file_half;

  reconstruction(int n_iter, std::string matrix, std::string mean):
    n_iter(n_iter) {
    ibstream in(matrix);

    if(!in) {
      std::cerr<<"error opening matrix file '"<<matrix<<"'"<<std::endl;
      exit(-1);
    }

    typedef uint32_t file_int;
    typedef uint16_t file_half;

    file_int in_magic;

    in >> in_magic;

    if (in_magic != detector_ring<F>::magic_f) {
      throw("invalid input system matrix file");
    }

    in >> n_pixels;
    in >> emissions;
    in >> n_tubes;

    total_n_pixels=n_pixels*n_pixels;
    rho_.resize(total_n_pixels,(F)0.0);
    rho_detected_.resize(total_n_pixels,(F)0.0);

    scale.resize(total_n_pixels,(F)0.0);

    std::ifstream mean_file(mean);
    if(!mean_file) {
      std::cerr<<"error opening mean file '"<<mean<<"'"<<std::endl;
      exit(-1);
    }
    std::vector<mean_per_lor> lor_mean;

    clock_t start=clock();
    for (;;) {

      int x,y,value;

      mean_file >> y >> x >> value;
      if (mean_file.eof()) break;

      mean_per_lor temp_obj;

      lor_type lor(x,y);

      temp_obj.lor = lor;
      temp_obj.n = value;

      lor_mean.push_back(temp_obj);
    }
    clock_t stop=clock();

    double time=static_cast<double>(stop-start)/CLOCKS_PER_SEC;
    std::cout<<"means read time = "<<time<<"s\n";

    std::vector<hits_per_pixel> pixels;

    int index = 0;
    n_non_zero_elements_=0;

    start=clock();
    for (;;) {

      file_half a, b;
      in >> a >> b;
      if (in.eof()) break;

      lor_type lor(a,b);

      n.push_back(get_mean_per_lor(a,b,lor_mean));

      file_int count;

      in >> count;

      for (int i = 0;i < count; ++i) {

        file_half x, y;
        file_int hits;

        in >> x >> y >> hits;
        hits_per_pixel data;
        data.probability = static_cast<F>(hits/static_cast<F>(emissions));
        data.index = LOCATION(x,y,n_pixels);
        scale[data.index] += data.probability;
        n_non_zero_elements_++;
        pixels.push_back(data);
      }

      system_matrix.push_back(pixels);
      pixels.clear();
      index++;
    }
    stop=clock();

    time=static_cast<double>(stop-start)/CLOCKS_PER_SEC;
    std::cout<<"matrix read time = "<<time<<"s\n";


    for (auto it_vector = system_matrix.begin();
         it_vector != system_matrix.end();
         it_vector++) {
      for (auto it_list = it_vector->begin();
           it_list != it_vector->end();
           it_list++) {
        int pixel=it_list->index;
        it_list->probability/=scale[pixel];
      }
    }

    
      for (int p = 0; p < n_pixels * n_pixels;++p) {
        if(scale[p]>0)
          rho_detected_[p]=(F)1.0;
      }
    std::cout
      << "   Pixels: " << n_pixels  << std::endl
      << "Emissions: " << emissions << std::endl
      << "Detectors: " << n_tubes   << std::endl
      << "     LORs: " << system_matrix.size() << std::endl;
    std::cout<< "Non zero elements: "<<n_non_zero_elements_<<std::endl;
  }

  ~reconstruction() {
  }


  void emt(int n_iter) {


    F y[n_pixels * n_pixels];
    F u;

    clock_t start=clock();

    for (int i = 0; i < n_iter; ++i) {
      std::cout << ".";
      std::cout.flush();

      for (int p = 0; p < n_pixels * n_pixels;++p) {
        y[p] = (F)0.0;
      }

      int t = 0;
      for (auto it_vector = system_matrix.begin();
           it_vector != system_matrix.end();
           it_vector++) {
        u = (F)0.0;
        if (n[t]> 0) {
          for (auto it_list = it_vector->begin();
               it_list != it_vector->end();
               it_list++) {
            u += rho_detected_[it_list->index] * it_list->probability;
          }
          F phi = n[t]/u;
          for (auto it_list = it_vector->begin();
               it_list != it_vector->end();
               ++it_list) {
            y[it_list->index] += phi * it_list->probability;
          }
        }
        t++;
      }

      for (int p = 0; p < n_pixels * n_pixels; ++p) {
        if (scale[p] > 0) {
          rho_detected_[p] *= y[p] ;
         
        }
      }
    }
    clock_t stop=clock();
    std::cout << std::endl;

    for (int p = 0; p < n_pixels * n_pixels; ++p) {
      if (scale[p] > 0) {
        rho_[p] = (rho_detected_[p]/scale[p]) ;
      }
    }

    double time=static_cast<double>(stop-start)/CLOCKS_PER_SEC;
    std::cout<<"time = "<<time<<"s time/iter = "<<time/n_iter<<"s\n";
    std::cout<<"op/sec = "<<n_non_zero_elements_*(n_iter/time)<<"\n";
  }

  int get_n_pixels() { return n_pixels; }

  F rho(int p) const {return rho_[p];}
  F rho_detected(int p) const {return rho_detected_[p];}
  std::vector<F> rho() const {return rho_;}
  std::vector<F> rho_detected() const {return rho_detected_;}

private:

  int n_tubes;
  int n_pixels;
  int total_n_pixels;
  int n_iter;
  int emissions;
  int n_lors;
  int n_non_zero_elements_;
  std::vector< std::vector<hits_per_pixel> > system_matrix;
  std::vector<int> list_of_lors;
  std::vector<int> n;
  std::vector<F>  scale;

  std::vector<F> rho_;
  std::vector<F> rho_detected_;

  int get_mean_per_lor(file_half &a,file_half &b,std::vector<mean_per_lor> &mean) {
    for (auto it = mean.begin(); it != mean.end(); ++it) {
      if (a == it->lor.first && b == it->lor.second) {
#if DEBUG
        std::cout << "Equal: lor(" << it->lor.first << "," << it->lor.second << ")" << std::endl;
#endif
        return it->n;
      }
    }

    return 0;
  }

};
