#pragma once

template <typename F> class ConstantErrorsModel {

 public:
  ConstantErrorsModel(F sigma_z, F sigma_l, F correlation = 0.0)
      : sigma_z_(sigma_z), sigma_l_(sigma_l), correlation_(correlation) {}

 private:
  F sigma_z_;
  F sigma_l_;
  F correlation_;
};
