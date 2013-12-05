#pragma once

#include <cstddef>

template <std::size_t MaxSize, typename T> class Array {
 public:
  Array() : s(0) {}
  template <typename... Ts> Array(Ts&&... e) : s(sizeof...(e)), v{ e... } {}

  // minimal std::vector compatibility
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef pointer iterator;
  typedef const_pointer const_iterator;
  typedef ptrdiff_t difference;
  typedef std::size_t size_type;

  iterator begin() { return v; }
  const_iterator begin() const { return v; }

  iterator end() { return v + s; }
  const_iterator end() const { return v + s; }

  size_type size() const { return s; }
  static size_type max_size() { return MaxSize; }

  void push_back(const value_type& val) { v[s++] = val; }

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
  T v[MaxSize];
};
