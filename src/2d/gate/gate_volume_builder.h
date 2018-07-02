#ifndef GATE_VOLUME_BUILDER_H
#define GATE_VOLUME_BUILDER_H

#include "gate_volume.h"

namespace Gate
{
namespace D2
{

template <typename F> Volume<F>* build_big_barrel_volume()
{
  using Box = Box<F>;
  using Vector = typename Box::Vector;
  using Cylinder = Cylinder<F>;

  Volume<F>* world = new Box(2, 2);

  auto layer_1 = new Cylinder(0.425 - 0.005, 0.425 + 0.005);
  world->attach_daughter(layer_1);
  auto scintillator_1 = new Box(0.021, 0.009);
  scintillator_1->set_translation(Vector(0.425, 0));
  scintillator_1->attach_repeater(new Gate::D2::Ring<F>(48, Vector(0, 0)));
  scintillator_1->attach_crystal_sd();
  layer_1->attach_daughter(scintillator_1);

  auto layer_2 = new Cylinder(0.4675 - 0.005, 0.4675 + 0.005);
  world->attach_daughter(layer_2);
  auto scintillator_2 = new Box(0.021, 0.009);
  scintillator_2->set_translation(Vector(0.4675, 0));
  scintillator_2->attach_repeater(
    new Gate::D2::Ring<F>(48, Vector(0, 0), M_PI / 48));
  scintillator_2->attach_crystal_sd();
  layer_2->attach_daughter(scintillator_2);

  auto layer_3 = new Cylinder(0.575 - 0.005, 0.575 + 0.005);
  world->attach_daughter(layer_3);
  auto scintillator_3 = new Box(0.021, 0.009);
  scintillator_3->set_translation(Vector(0.575, 0));
  scintillator_3->attach_repeater(
    new Gate::D2::Ring<F>(96, Vector(0, 0), M_PI / 96));
  scintillator_3->attach_crystal_sd();
  layer_3->attach_daughter(scintillator_3);

  return world;
}

template <typename F> Volume<F>* build_new_full_scanner_volume()
{
  using Box = Gate::D2::Box<F>;
  using Vector = typename Box::Vector;
  using Cylinder = Gate::D2::Cylinder<F>;

  auto world = new Box(2, 2);

  auto layer_new = new Cylinder(0.35, 0.4);
  world->attach_daughter(layer_new);

  auto module = new Box(0.026, 0.0085);
  module->set_translation(Vector(0.37236, 0));
  module->attach_repeater(new Gate::D2::Ring<F>(24, Vector(0.0, 0.0)));

  layer_new->attach_daughter(module);

  auto scintillator = new Box(0.024, 0.006);
  scintillator->attach_repeater(new Gate::D2::Linear<F>(13, Vector(0, 0.007)));
  scintillator->attach_crystal_sd();

  module->attach_daughter(scintillator);

  auto layer_1 = new Cylinder(0.425 - 0.005, 0.425 + 0.005);
  world->attach_daughter(layer_1);
  auto scintillator_1 = new Box(0.021, 0.009);
  scintillator_1->set_translation(Vector(0.425, 0));
  scintillator_1->attach_repeater(new Gate::D2::Ring<F>(48, Vector(0, 0)));
  scintillator_1->attach_crystal_sd();
  layer_1->attach_daughter(scintillator_1);

  auto layer_2 = new Cylinder(0.4675 - 0.005, 0.4675 + 0.005);
  world->attach_daughter(layer_2);
  auto scintillator_2 = new Box(0.021, 0.009);
  scintillator_2->set_translation(Vector(0.4675, 0));
  scintillator_2->attach_repeater(
    new Gate::D2::Ring<F>(48, Vector(0, 0), M_PI / 48));
  scintillator_2->attach_crystal_sd();
  layer_2->attach_daughter(scintillator_2);

  auto layer_3 = new Cylinder(0.575 - 0.005, 0.575 + 0.005);
  world->attach_daughter(layer_3);
  auto scintillator_3 = new Box(0.021, 0.009);
  scintillator_3->set_translation(Vector(0.575, 0));
  scintillator_3->attach_repeater(
    new Gate::D2::Ring<F>(96, Vector(0, 0), M_PI / 96));
  scintillator_3->attach_crystal_sd();
  layer_3->attach_daughter(scintillator_3);

  return world;
}

template <typename F> Volume<F>* build_new_module_scanner_volume()
{
  using Box = Gate::D2::Box<F>;
  using Vector = typename Box::Vector;
  using Cylinder = Gate::D2::Cylinder<F>;

  auto world = new Box(2, 2);

  auto layer_new = new Cylinder(0.35, 0.4);
  world->attach_daughter(layer_new);

  auto module = new Box(0.026, 0.0085);
  module->set_translation(Vector(0.37236, 0));
  module->attach_repeater(new Gate::D2::Ring<F>(24, Vector(0.0, 0.0)));

  layer_new->attach_daughter(module);

  auto scintillator = new Box(0.024, 0.006);
  scintillator->attach_repeater(new Gate::D2::Linear<F>(13, Vector(0, 0.007)));
  scintillator->attach_crystal_sd();

  module->attach_daughter(scintillator);

  return world;
}
}  // namespace D2
}  // namespace Gate

#endif  // GATE_VOLUME_BUILDER_H
