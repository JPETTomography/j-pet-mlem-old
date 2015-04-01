#pragma once
#include <initializer_list>
#include <valarray>

const double Degree = M_PI / 180.0;

namespace Geometry {

  namespace Base {

    template <int N, typename F> class Point;

    template <int N, typename F = double> class Vector {
     public:
      Vector() {
        for (int i = 0; i < N; i++)
          v_[i] = F();
      }

      Vector(F v[N]) {
        for (int i = 0; i < N; i++)
          v_[i] = v[i];
      }

      Vector(std::initializer_list<F> v) {
        auto it = v.begin();
        for (int i = 0; i < N; i++)
          v_[i] = it[i];
      }

      Vector& operator+=(const Vector& rhs) {
        for (int i = 0; i < N; i++)
          v_[i] += rhs.v_[i];
        return *this;
      }

      Vector& operator-=(const Vector& rhs) {
        for (int i = 0; i < N; i++)
          v_[i] -= rhs.v_[i];
        return *this;
      }

      Vector& operator+=(F rhs) {
        for (int i = 0; i < N; i++)
          v_[i] += rhs;
        return *this;
      }

      Vector& operator-=(F rhs) {
        for (int i = 0; i < N; i++)
          v_[i] -= rhs;
        return *this;
      }

      Vector& operator*=(F rhs) {
        for (int i = 0; i < N; i++)
          v_[i] *= rhs;
        return *this;
      }

      Vector& operator/=(F rhs) {
        for (int i = 0; i < N; i++)
          v_[i] /= rhs;
        return *this;
      }

      F& operator[](int i) { return v_[i]; }
      F operator[](int i) const { return v_[i]; }

      friend class Point<N, F>;

     protected:
      F v_[N];
    };

    template <int N, typename F> class Point {
     public:
      Point() {
        for (int i = 0; i < N; i++)
          p_[i] = F();
      }
      Point(F p[N]) {
        for (int i = 0; i < N; i++)
          p_[i] = p[i];
      }

      Point(std::initializer_list<F> p) {
        int i = 0;
        auto it = p.begin();
        for (; i < N; ++i, ++it)
          p_[i] = *it;
      }

      Point& operator+=(const Vector<N, F>& rhs) {
        for (int i = 0; i < N; i++)
          p_[i] += rhs.v_[i];
        return *this;
      };
      Point& operator-=(const Vector<N, F>& rhs) {
        for (int i = 0; i < N; i++)
          p_[i] -= rhs.v_[i];
        return *this;
      }

      Point& operator+=(F rhs) {
        for (int i = 0; i < N; i++)
          p_[i] += rhs;
        return *this;
      };

      Point& operator-=(F rhs) {
        for (int i = 0; i < N; i++)
          p_[i] -= rhs;
        return *this;
      }

      Point& operator*=(F rhs) {
        for (int i = 0; i < N; i++)
          p_[i] *= rhs;
        return *this;
      };

      Point& operator/=(F rhs) {
        for (int i = 0; i < N; i++)
          p_[i] /= rhs;
        return *this;
      }

      //      Vector<N,F> operator-(const Point& rhs) {
      //        return Vector<N,F>(this->p_-rhs.p_);
      //      }

      F& operator[](int i) { return p_[i]; }
      F operator[](int i) const { return p_[i]; }

     protected:
      F p_[N];
    };
  }

  template <int N, typename F = double>
  class Vector : public Base::Vector<N, F> {
   public:
    Vector(F v[N]) : Base::Vector<N, F>(v) {}
  };
  template <int N, typename F = double> class Point : public Base::Point<N, F> {
   public:
    Point(F v[N]) : Base::Point<N, F>(v) {}
  };

  template <typename F> class Vector<2, F> : public Base::Vector<2, F> {
   public:
    Vector() : Base::Vector<2, F>(), x(this->v_[0]), y(this->v_[1]) {}
    Vector(std::initializer_list<F> v)
        : Base::Vector<2, F>(v), x(this->v_[0]), y(this->v_[1]) {}
    Vector(F x_a, F y_a)
        : Base::Vector<2, F>(), x(this->v_[0]), y(this->v_[1]) {
      this->v_[0] = x_a;
      this->v_[1] = y_a;
    }

    Vector(F v[2]) : Base::Vector<2, F>(v), x(this->v_[0]), y(this->v_[1]) {}

    Vector(const Vector& other) : x(this->v_[0]), y(this->v_[1]) {
      for (int i = 0; i < 2; ++i)
        this->v_[i] = other[i];
    }

    Vector& rotate(F angle) {
      F s = sin(angle);
      F c = cos(angle);
      F x = (*this)[0];
      F y = (*this)[1];

      (*this)[0] = x * c - y * s;
      (*this)[1] = y * c + x * s;
    }

    Vector& operator=(const Vector& rhs) {
      for (int i = 0; i < 2; i++)
        this->v_[i] = rhs.v_[i];
      return *this;
    }

    F& x;
    F& y;
  };

  template <int N, typename F = double>
  Vector<N, F> operator+(const Vector<N, F>& v1, const Vector<N, F>& v2) {
    Vector<N, F> r(v1);
    r += v2;
    return r;
  }

  template <int N, typename F = double>
  Vector<N, F> operator-(const Vector<N, F>& v1, const Vector<N, F>& v2) {
    Vector<N, F> r(v1);
    r -= v2;
    return r;
  }

  template <int N, typename F = double>
  Vector<N, F> operator*(const Vector<N, F>& v1, F s) {
    Vector<N, F> r(v1);
    r *= s;
    return r;
  }

  template <int N, typename F = double>
  Vector<N, F> operator*(F s, const Vector<N, F>& v1) {
    Vector<N, F> r(v1);
    r *= s;
    return r;
  }

  template <int N, typename F = double>
  Vector<N, F> operator/(const Vector<N, F>& v1, F s) {
    Vector<N, F> r(v1);
    r /= s;
    return r;
  }

  template <typename F> class Point<2, F> : public Base::Point<2, F> {
   public:
    Point() : Base::Point<2, F>(), x(this->p_[0]), y(this->p_[1]) {}
    Point(F v[2]) : Base::Point<2, F>(v), x(this->p_[0]), y(this->p_[1]) {}
    Point(std::initializer_list<F> p)
        : Base::Point<2, F>(p), x(this->p_[0]), y(this->p_[1]) {}
    Point(F x_a, F y_a) : Base::Point<2, F>(), x(this->p_[0]), y(this->p_[1]) {
      this->p_[0] = x_a;
      this->p_[1] = y_a;
    }

    F& x;
    F& y;

    Point(const Point& other) : x(this->p_[0]), y(this->p_[1]) {
      for (int i = 0; i < 2; i++)
        this->p_[i] = other[i];
    }

    Point& operator=(const Point& rhs) {
      for (int i = 0; i < 2; i++)
        this->p_[i] = rhs.p_[i];
      return *this;
    }

    Point& rotate(F angle, Point center) {
      Vector<2, F> v = (*this) - center;
      v.rotate(angle);
      Point p = center + v;
      (*this)[0] = p[0];
      (*this)[1] = p[1];

      return *this;
    }
    Point& rotate(F angle) {
      Vector<2, F> v;
      v[0] = (*this)[0];
      v[1] = (*this)[1];
      v.rotate(angle);
      (*this)[0] = v[0];
      (*this)[1] = v[1];

      return *this;
    }
  };

  template <int N, typename F = double>
  Point<N, F> operator+(const Vector<N, F>& v, const Point<N, F>& p) {
    Point<N, F> r(p);
    r += v;
    return r;
  }

  template <int N, typename F = double>
  Point<N, F> operator+(const Point<N, F>& p, const Vector<N, F>& v) {
    Point<N, F> r(p);
    r += v;
    return r;
  }

  template <int N, typename F = double>
  Point<N, F> operator-(const Point<N, F>& p, const Vector<N, F>& v) {
    Point<N, F> r(p);
    r -= v;
    return r;
  }

  template <int N, typename F = double>
  Vector<N, F> operator-(const Point<N, F>& p1, const Point<N, F>& p2) {
    Vector<N, F> v;
    for (int i = 0; i < N; ++i)
      v[i] = p1[i] - p2[i];

    return v;
  }

  template <int N, typename F = double>
  Vector<N, F> operator+(const Vector<N, F>& v1, F s);
  template <int N, typename F = double>
  Vector<N, F> operator+(F s, const Vector<N, F>& v1);
}
