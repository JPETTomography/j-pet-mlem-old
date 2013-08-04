#ifndef SCINTILLATOR_H
#define SCINTILLATOR_H

#endif  // SCINTILLATOR_H

template <typename T = float>
class scintillator {

 private:
  T y;
  T z;
  T l;

 public:
  scintillator(T y, T z, T l) : y(y), z(z), l(l) {}
  T get_y() const { return y; }
  T get_z() const { return z; }
  T get_l() const { return l; }
};
