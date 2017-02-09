#ifndef GATE_VOLUME_H
#define GATE_VOLUME_H
namespace Gate {

class Repeater {};

class Ring : public Repeater {};

class Linear : public Repeater {};

class Volume {
 public:
  virtual void build() = 0;
};

class Box : public Volume {};

class Cylinder : public Volume {};
}

#endif  // GATE_VOLUME_H
