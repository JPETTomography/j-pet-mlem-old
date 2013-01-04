#pragma once

template <typename F = double, typename HitType = int>
class TriangularPixelMap {
public:
  typedef HitType hit_type;
  typedef hit_type *pixels_type;
  typedef uint8_t bitmap_pixel_type;
  typedef typename detector_ring<F>::lor_type lor_type;

      // reserve for pixel stats
  TriangularPixelMap(int n_pixels_a)
 :n_pixels(n_pixels_a),
  n_pixels_2(n_pixels/2),
  n_t_matrix_pixels( n_pixels_a/2 * (n_pixels_a/2+1) / 2 ) {
    t_hits = new hit_type[n_t_matrix_pixels]();
  }

  ~TriangularPixelMap() {
    delete [] t_hits;
  }

  static constexpr size_t t_pixel_index(size_t x, size_t y) {
    return y*(y+1)/2 + x;
  }

  /**Computes pixel index and determines symmetry number based on pixel
     position
     @param x        pixel x coordinate (0..n_pixels)
     @param x        pixel y coordinate (0..n_pixels)
     @param diag     outputs true if abs(x)==abs(y)
     @param symmetry outputs symmetry number (0..7)
  */
  size_t pixel_index(ssize_t x, ssize_t y, bool &diag, int &symmetry) const {
    // shift so 0,0 is now center
    x -= n_pixels_2; y -= n_pixels_2;
    // mirror
    symmetry = 0;
    if (x < 0) {
      x = -x-1;
      symmetry |= 2;
    };
    if (y < 0) {
      y = -y-1;
      symmetry |= 1;
    }
    // triangulate
    if (x > y) {
      std::swap(x, y);
      symmetry |= 4;
    };
    diag = (x == y);
    return t_pixel_index(x, y);
  }

  hit_type hits(size_t x, size_t y) const {
    bool diag; int symmetry;
    auto i_pixel = pixel_index(x, y, diag, symmetry);
    return this->t_hits[i_pixel] * (diag ? 2 : 1);
  }


  template<class FileWriter>
  void output_bitmap(FileWriter &fw) {
    fw.template write_header<bitmap_pixel_type>(n_pixels, n_pixels);
    hit_type pixel_max = 0;
    for (auto y = 0; y < n_pixels; ++y) {
      for (auto x = 0; x < n_pixels; ++x) {
        pixel_max = std::max(pixel_max, hits(x, y));
      }
    }
    auto gain = static_cast<double>(std::numeric_limits<bitmap_pixel_type>::max()) / pixel_max;
    for (auto y = 0; y < n_pixels; ++y) {
      bitmap_pixel_type row[n_pixels];
      for (auto x = 0; x < n_pixels; ++x) {
        row[x] = std::numeric_limits<bitmap_pixel_type>::max() - gain * hits(x, y);
      }
      fw.write_row(row);
    }
  }



protected:
  pixels_type t_hits;
  int n_pixels; 
  int n_pixels_2; 
  int n_t_matrix_pixels;
private:
};
