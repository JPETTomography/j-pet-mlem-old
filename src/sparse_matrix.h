#pragma once

template <typename LORType, typename SType = int, typename HitType = int>
class SparseMatrix : public std::vector<std::tuple<LORType, SType, HitType>> {
  typedef std::vector<std::tuple<LORType, SType, HitType>> Super;

 public:
  typedef LORType LOR;
  typedef SType S;
  typedef HitType Hit;
  typedef std::tuple<LOR, S, Hit> Element;

  void sort_by_lor() { std::sort(Super::begin(), Super::end(), SortByLOR()); }
  void sort_by_pixel() {
    std::sort(Super::begin(), Super::end(), SortByPixel());
  }

 private:
  struct SortByPixel {
    bool operator()(const Element& a, const Element& b) const {
      return std::get<1>(a) < std::get<1>(b);
    }
  };

  struct SortByLOR {
    bool operator()(const Element& a, const Element& b) const {
      return std::get<0>(a) < std::get<0>(b);
    }
  };
};
