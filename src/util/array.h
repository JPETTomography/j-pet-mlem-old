#pragma once

#include <new>
#include <type_traits>
#include <utility>

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
  array() : s(0) {}

#if !_MSC_VER || __CUDACC__
  template <typename... Args>
  array(Args&&... e)
      : s(sizeof...(e)), v{ *reinterpret_cast<storage_type*>(&e)... } {}
#else
 private:
  /// Emulates initialization using recursive variadic templates
  template <std::size_t I, typename Arg, typename... Args>
  void _init(Arg&& first, Args&&... rest) {
    new (&v[I]) value_type(std::forward<Arg&&>(first));
    _init<I + 1>(std::forward<Args&&>(rest)...);
  }
  template <std::size_t I, typename Arg> void _init(Arg&& last) {
    new (&v[I]) value_type(std::forward<Arg&&>(last));
  }

 public:
  template <typename... Args> array(Args&&... e) : s(sizeof...(e)) {
    _init<0>(std::forward<Args&&>(e)...);
  }
#endif

  /// Returns if the array is full (has max number of elements)
  bool full() const { return (s == MaxSize); }

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

  iterator begin() { return reinterpret_cast<pointer>(v); }
  const_iterator begin() const { return reinterpret_cast<const_pointer>(v); }

  iterator end() { return reinterpret_cast<pointer>(v + s); }
  const_iterator end() const { return reinterpret_cast<const_pointer>(v + s); }

  size_type size() const { return s; }
  static size_type max_size() { return MaxSize; }

  void push_back(const value_type& val) { new (&v[s++]) value_type(val); }

  template <typename... Args> void emplace_back(Args&&... args) {
    new (&v[s++]) value_type(std::forward<Args&&>(args)...);
  }

  reference at(size_type i) { return *reinterpret_cast<pointer>(&v[i]); }
  const_reference at(size_type i) const {
    return *reinterpret_cast<const_pointer>(&v[i]);
  }

  reference operator[](size_type i) { return at(i); }
  const_reference operator[](size_type i) const { at(i); }

  reference front() { return at(0); }
  const_reference front() const { at(0); }

  reference back() { return at(s - 1); }
  const_reference back() const { return at(s - 1); }

 private:
  std::size_t s;
  storage_type v[MaxSize];
};
}  // util
