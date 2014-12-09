#pragma once

#include <new>
#include <type_traits>
#include <utility>

#ifndef _
#include "cuda/compat.h"
#endif

namespace util {

/// Stack based replacement for \c std::vector

/// This class is drop-in replacement for \c std::vector for all cases when
/// maximum size is known at compile time.
/// This way we can pass everything via stack, omitting unnecessary allocations.
/// It uses internal storage type and user specified alignment.
template <std::size_t MaxSize,                 //< maximum size
          typename T,                          //< type carried by
          std::size_t Alignment = alignof(T),  //< storage alignment
          typename S =                         //< storage type
          /// (automatically determined only if T trivially_destructible)
          typename std::enable_if<
              std::is_trivially_destructible<T>::value,
              std::aligned_storage<sizeof(T), Alignment>>::type::type>
class array {
  typedef S storage_type;  ///< must be same size and alignment as value_type

 public:
  /// Make an empty array.
  _ array() : s(0) {}
  /// Make one element array with given value.
  explicit _ array(T&& arg) : s(1) {
    new (&v[0]) value_type(std::forward<T>(arg));
  }
  /// Make n-element array with given values.
  template <typename... Args>
  explicit _ array(T&& arg, Args&&... args)
      : s(sizeof...(args) + 1) {
    __init<0>(std::forward<T>(arg), std::forward<Args>(args)...);
  }

  /// Returns if the array is full (has max number of elements)
  _ bool full() const { return (s == MaxSize); }

  // minimal std::vector compatibility
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef std::size_t size_type;
  typedef pointer iterator;
  typedef const_pointer const_iterator;

  enum : std::size_t {
    value_size = sizeof(value_type),
    storage_size = sizeof(storage_type),
    alignment = Alignment
  };

  _ iterator begin() { return reinterpret_cast<pointer>(v); }
  _ const_iterator begin() const { return reinterpret_cast<const_pointer>(v); }

  _ iterator end() { return reinterpret_cast<pointer>(v + s); }
  _ const_iterator end() const {
    return reinterpret_cast<const_pointer>(v + s);
  }

  _ size_type size() const { return s; }
  _ static size_type max_size() { return MaxSize; }

  _ void push_back(const value_type& val) { new (&v[s++]) value_type(val); }

  template <typename... Args> _ void emplace_back(Args&&... args) {
    new (&v[s++]) value_type(std::forward<Args&&>(args)...);
  }

  _ reference at(size_type i) { return reinterpret_cast<reference>(v[i]); }
  _ const_reference at(size_type i) const {
    return reinterpret_cast<const_reference>(v[i]);
  }

  _ reference operator[](size_type i) { return at(i); }
  _ const_reference operator[](size_type i) const { return at(i); }

  _ reference front() { return at(0); }
  _ const_reference front() const { return at(0); }

  _ reference back() { return at(s - 1); }
  _ const_reference back() const { return at(s - 1); }

  _ bool operator==(const array& other) const {
    if (other.size() != s)
      return false;
    for (size_type i = 0; i < s; ++i) {
      if (other[i] != at(i))
        return false;
    }
    return true;
  }

 private:
  storage_type v[MaxSize];
  std::size_t s;

  template <std::size_t I, typename... Args>
  _ void __init(T&& first, Args&&... rest) {
    new (&v[I]) value_type(std::forward<T>(first));
    __init<I + 1>(std::forward<Args>(rest)...);
  }
  template <std::size_t I> _ void __init(T&& last) {
    new (&v[I]) value_type(std::forward<T>(last));
  }
};
}  // util

#ifdef TEST_CASE
namespace Catch {
template <std::size_t MaxSize, typename T, std::size_t Alignment, typename S>
struct StringMaker<util::array<MaxSize, T, Alignment, S>> {
  using Array = util::array<MaxSize, T, Alignment, S>;
  using size_type = typename Array::size_type;
  static std::string convert(const Array& array) {
    std::ostringstream oss;
    oss << "[";
    for (size_type i = 0; i < array.size(); ++i) {
      if (i)
        oss << ", ";
      oss << array[i];
    }
    oss << "]";
    return oss.str();
  }
};
}
#endif
