#pragma once


/**
   This class represents a system matrix that stores the content in
   "pixel major" mode. That is for each pixel a list w lors is kept.

   The idea behind it is that the MC simulations are done on per pixel basis.
   That means that in the same time only few pixels are processed in parallel.
   On shiva for example the biggest machine has 24 cores.
   We can afford to alloc memory for lors in an uneficient way providing for
   quick acces and reduce them after the simulations of  this pixel is finished.

   It alos means that as different threads are processing different pixels,
   there is no need for synchronisation in add_pixels.

   The reconstruction is however done using the lor_major matrix.
   So this class has the possibility to write the matrix down in the
   triangular lor_major and full lor_major format.

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
     n_entries_(0),
     index_to_lor_(SimpleLor::n_lors()) {
    for(auto l_it=SimpleLor::begin();l_it!=SimpleLor::end();++l_it) {
      auto lor=*l_it;
      index_to_lor_[t_lor_index(lor)]=lor;
    }

  };


  int n_entries() const {return n_entries_;}
  int n_lors(int p) const {return pixel_count_[p];}
  int total_n_pixels() const {return total_n_pixels_;}

  void add_to_t_matrix(const LorType &lor, int i_pixel) {

    if(!pixel_tmp_[i_pixel]) {
      pixel_tmp_[i_pixel]=new int[n_lors_]();
    }
    if(pixel_tmp_[i_pixel][t_lor_index(lor)]==0){
      pixel_count_[i_pixel]++;
      n_entries_++;
    }
    pixel_tmp_[i_pixel][t_lor_index(lor)]++;
  };

  ~PixelMajorSystemMatrix() {
    for(int p=0;p<total_n_pixels_;++p) {
      if(pixel_tmp_[p]) {
        delete [] pixel_tmp_[p];
      }
    }
  }

  int t_get_element(LorType lor, int i_pixel) {
    auto hit=std::lower_bound(pixel_[i_pixel].begin(),
                              pixel_[i_pixel].end(),
                              std::make_pair(lor,0),
                              HitTypeComparator());

    if(hit== pixel_[i_pixel].end())
      return 0;
    return (*hit).second;
  }

  void finalize_pixel(int i_pixel) {
    int count=0;
    pixel_[i_pixel].resize(pixel_count_[i_pixel]);
    int it=0;
    for(int lor=0;lor<n_lors_;++lor) {
      int hits;
      if((hits=pixel_tmp_[i_pixel][lor])>0) {
        pixel_[i_pixel][it]=std::make_pair(index_to_lor_[lor],hits);
        it++;
      }
    }
    delete [] pixel_tmp_[i_pixel];
    pixel_tmp_[i_pixel]=(int *)0;
    std::sort(pixel_[i_pixel].begin(),
              pixel_[i_pixel].end(),
              HitTypeComparator()
              );
  }

  static size_t t_lor_index(const LorType &lor) {
    return LorType::t_index(lor);
  }

  typedef std::pair< std::pair<LorType,int>,int> PairType;

  void to_pairs() {
    pair_.reserve(n_entries());
    for(int p=0;p<total_n_pixels();++p) {
      for(auto it=pixel_[p].begin();it!=pixel_[p].end();++it) {
        pair_.push_back(
                        std::make_pair(std::make_pair( (*it).first,p ),
                                       (*it).second)
                        );
      }
    }
  }

  PairType pair(int p) const {return pair_[p];} 

  void sort_pairs_by_lors() {
    std::sort(pair_.begin(),pair_.end(),LorSorter());
  };
  void sort_pairs_by_pixels() {
    std::sort(pair_.begin(),pair_.end(),PixelSorter());
  };


private:

  struct HitTypeComparator {
    bool operator()(const HitType &a,const HitType &b) const {
      return LorType::less(a.first,b.first);
    }
  };

  struct LorSorter {
    bool operator()(const PairType &a,const PairType &b) const {
      return LorType::less(a.first.first,b.first.first);
    }
  };

  struct PixelSorter {
    bool operator()(const PairType &a,const PairType &b) const {
      return a.first.second<b.first.second;
    }
  };

  int n_pixels_;
  int n_pixels_half_;
  int total_n_pixels_;
  int n_lors_;
  int n_entries_;
  std::vector<int*> pixel_tmp_;
  std::vector< std::vector<HitType> > pixel_;
  std::vector<int> pixel_count_;
  std::vector<LorType> index_to_lor_;
  std::vector< PairType > pair_;

};
