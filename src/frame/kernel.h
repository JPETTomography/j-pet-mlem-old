#pragma once

template <typename F,
          template <typename FT> class Detector,
          template <typename FT> class ErrorsModel>
class Kernel {
 public:
  typedef Detector<F> DetectorType;
  typedef ErrorsModel<F> ErrorsModelType;
  Kernel(Detector<F> detector, ErrorsModel<F> errors_model)
      : detector_(detector), errors_model_(errors_model) {}

  class EventKernel {
   public:
    /// NYI
    F operator()(F dy __attribute__((unused)), F dz __attribute__((unused))) {
      return 0.0;
    }
  };

  typedef EventKernel EventKernelType;

  /// NYI
  EventKernel MakeEventKernel(EventImageAngle<F> event
                              __attribute__((unused))) {
    return EventKernel();
  }

 private:
  DetectorType detector_;
  ErrorsModelType errors_model_;
};
