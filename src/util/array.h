#pragma once

#include <cstddef>
#if _MSC_VER
#include <array>
#endif

namespace util {

/// Stack based replacement for \c std::vector

/// This class is drop-in replacement for \c std::vector for all cases when
/// maximum size is known at compile time.
/// This way we can pass everything via stack, omitting unnecessary allocations.
template <std::size_t MaxSize, typename T> class array {
 public:
#if !_MSC_VER
  array() : s(0), v() {}
#endif
#if __INTEL_COMPILER
  template <typename... Ts> array(Ts... e) : s(sizeof...(e)), v{ e... } {}
#elif _MSC_VER
  template <typename... Ts>
  array(Ts... e)
      : s(sizeof...(e)), v(std::array<T, MaxSize>{ e... }) {}
#else
  template <typename... Ts> array(Ts&&... e) : s(sizeof...(e)), v{ e... } {}
#endif

  /// returns if the array is full (has max number of elements)
  bool full() const { return (s == MaxSize); }

  // minimal std::vector compatibility
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef ptrdiff_t difference;
  typedef std::size_t size_type;

#if !_MSC_VER || __CUDACC__
  typedef pointer iterator;
  typedef const_pointer const_iterator;

  iterator begin() { return v; }
  const_iterator begin() const { return v; }

  iterator end() { return v + s; }
  const_iterator end() const { return v + s; }
#else  // MSC_VER
  typedef typename std::array<T, MaxSize>::iterator iterator;
  typedef typename std::array<T, MaxSize>::const_iterator const_iterator;

  iterator begin() { return v.begin(); }
  const_iterator begin() const { return v.begin(); }

  iterator end() { return v.end(); }
  const_iterator end() const { return v.end(); }
#endif

  size_type size() const { return s; }
  static size_type max_size() { return MaxSize; }

  void push_back(const value_type& val) { v[s++] = val; }

  template <typename... Args> void emplace_back(Args&&... args) {
    new (&v[s++]) value_type(std::forward<Args>(args)...);
  }

  reference at(size_type i) { return v[i]; }
  const_reference at(size_type i) const { return v[i]; }

  reference operator[](size_type i) { return v[i]; }
  const_reference operator[](size_type i) const { return v[i]; }

  reference front() { return v[0]; }
  const_reference front() const { return v[0]; }

  reference back() { return v[s - 1]; }
  const_reference back() const { return v[s - 1]; }

 private:
  std::size_t s;
#if !_MSC_VER
  T v[MaxSize];
#else
  std::array<T, MaxSize> v;
#endif
};
}  // util
