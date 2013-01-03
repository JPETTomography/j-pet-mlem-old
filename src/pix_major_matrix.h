#pragma once


/**
   This class represents a system matrix that stores the content in
   "pixel major" mode. That is for each pixel a list w lors is kept.

   The idea behind it is that the MC simulations are done on per pixel basis.
   That means that in the same time only few pixels are processed in parallel.
   On shiva for example the biggest machine has 24 cores.
   We can afford to alloc memory for lors in an uneficient way providing for
   quick acces and reduce them after the simulations for this pixel is finished.
*/


template<typename LorType, typename F = double>
class PixelMajorSystemMatrix {

public:
  typedef std::pair<LorType,int> HitType;


  PixelMajorSystemMatrix(detector_ring<F> dr_a,
                         int n_pixels_a,
                         F pixel_size_a)
    :n_lors_(SimpleLor::n_lors()),
     n_pixels_(n_pixels_a),
     n_pixels_half_(n_pixels_/2),
     total_n_pixels_(n_pixels_half_*(n_pixels_half_+1)/2),
     pixel_tmp_(total_n_pixels_,(int *)0),
     pixel_(total_n_pixels_),
     pixel_count_(total_n_pixels_),
     index_to_lor_(SimpleLor::n_lors()) {};

  void add_to_t_matrix(const LorType &lor, int i_pixel) {

    std::cerr<<"adding to  pixel"<<std::endl;

    if(!pixel_tmp_[i_pixel]) {
      std::cerr<<"alcating memory for "<<i_pixel<<" "<<n_lors_<<std::endl;
      pixel_tmp_[i_pixel]=new int[n_lors_]();
    }

    pixel_tmp_[i_pixel][t_lor_index(lor)]++;
    pixel_count_[i_pixel]++;
  };

  ~PixelMajorSystemMatrix() {
    for(int p=0;p<total_n_pixels_;++p) {
      if(pixel_tmp_[p])
        delete [] pixel_tmp_[p];
    }
  }

  int t_get_element(LorType lor, int i_pixel) {
    auto hit=std::lower_bound(pixel_[i_pixel].begin(),
                              pixel_[i_pixel].end(),
                              make_pair(lor,0),
                              HitTypeComparator());

    if(hit== pixel_[i_pixel].end())
      return 0;
    return hit.second;
  }

  void finalize_pixel(int i_pixel) {
    int count=0;
    pixel_[i_pixel].resize(pixel_count_[i_pixel]);
    int it=0;
    for(int lor=0;lor<n_lors_;++lor) {
      int hits;
      if(hits=pixel_tmp_[i_pixel][lor]>0) {
        pixel_[i_pixel][it]=std::make_pair(index_to_lor_[lor],hits);
        it++;
      }
    }
    delete [] pixel_tmp_[i_pixel];

    std::sort(pixel_[i_pixel].begin(),
              pixel_[i_pixel].end(),
              LorType::less
              );
  }

  static size_t t_lor_index(const LorType &lor) {
    return LorType::t_index(lor);
  }

private:

  struct HitTypeComparator {
    bool operator()(const HitType &a,const HitType &b) {
      return LorType::less(a.first,b.first);
    }
  };

  int n_pixels_;
  int n_pixels_half_;
  int total_n_pixels_;
  int n_lors_;
  std::vector<int*> pixel_tmp_;
  std::vector< std::vector<HitType> > pixel_;
  std::vector<int> pixel_count_;
  std::vector<LorType> index_to_lor_;
};
