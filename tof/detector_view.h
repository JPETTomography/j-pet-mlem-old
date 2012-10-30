#pragma once

#include "tof_detector.h"

template<typename F>
class DetectorView {

public:
  typedef F float_t;
  typedef ToF_Detector_2D<F> detector_t;

  DetectorView(GeometryPlot *gp, detector_t *detector):
    gp_(gp), detector_(detector) {};

  void render() {
    gp_->renderZYRectangle(-250, -360, 250, -340, 0, false);
    gp_->renderZYRectangle(-250, 340, 250, 360, 0, false);
  }

private:
   GeometryPlot *gp_;
   detector_t *detector_;
};
