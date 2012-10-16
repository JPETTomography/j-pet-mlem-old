#pragma once

template <typename F>
void sincos(F a, F &s, F &c) {
  s = sin(a);
  c = cos(a);
};

template<typename F> class event {
public:
  event(F x, F y, F phi)
  : x_(x)
  , y_(y)
  , phi_(phi) {};

  F x()   const { return x_; }
  F y()   const { return y_; }
  F phi() const { return phi_; }

private:
  F x_, y_;
  F phi_;
};
