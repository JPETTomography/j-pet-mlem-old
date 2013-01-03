#pragma once


/**
   This class represents a system matrix that stores the content in
   "pixel major" mode. That is for each pixel a list w lors is kept.

   The idea behind it is that the MC simulations are done on per pixel basis.
   That means that in the same time only few pixels are processed in parallel.
   On shiva for example the biggest machine has 24 cores.
   We can afford to alloc memory for lors in an uneficient way providing for
   quick acces and reduce them after the simulations for this pixel is reduced.

*/


template<typename LorType>
class PixelMajorSystemMatrix {

public:
  typedef std::pair<LorType,int> HitType;

  void add_to_t_matrix(const LorType &lor, int i_pixel) {
    if(pixel_tmp_[i_pixel]==NULL) {
      pixel_tmp_[i_pixel]=new int[n_lors_]();
    }

    pixel_tmp_[i_pixel][t_lor_index(lor)]++;
    pixel_count_[i_pixel]++;
  };

  int t_get_element(lor_type lor, int i_pixel) {
    std::binary_search(pixel_[i_pixel].begin(),
                       pixel_[i_pixel].end(),
                       make_pair(lor,0),
                       HitTypeComparator());
  }

  void finalize_pixel(int i_pixel) {
    int count=0;
    pixel_[i_pixel].resize(pixel_counts_[i_pixel]);
    int it;
    for(int lor=0;i<n_lors_;++lor) {
      int hits;
      if(hits=pixel_tmp_[i_pixel][lor]>0) {
        pixel_[i_pixel][it]=std::make_pair<index_to_lor_[lor],hits>;
        it++;
      }
    }
    delete [] pixel[i_pixel];

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
    bool operator()(const HitType &a,const HitType &) {
      return LorType::less(a.first,b.first);
    }
  }
  int n_pixels_half_;
  int n_lors_;
  std::vector<*int> pixel_tmp_;
  std::vector< std::vector<HitType> > pixel_;
  std::vector<int> pixel_count_;
  std::vector<LorType> index_to_lor_;
};
