#pragma once

#include "tof_event.h"
#include "topet_simulator.h"

template <typename F> class EventView {
 public:
  typedef F float_t;
  typedef ToF_Event2D<F> event_t;
  typedef ToF_Track2D<F> track_t;
  typedef ToF_Detector_2D<F> detector_t;
  typedef typename TOPETSimulator<F>::const_iterator const_iterator;
  EventView(GeometryPlot* gp,
            detector_t* detector,
            const_iterator first,
            const_iterator last)
      : gp_(gp),
        detector_(detector),
        two_sets_(false),
        two_sets_mode_(false),
        first_(first),
        last_(last) {

  }

  EventView(GeometryPlot* gp,
            detector_t* detector,
            const_iterator first,
            const_iterator last,
            const_iterator first2)
      : gp_(gp),
        detector_(detector),
        two_sets_(true),
        two_sets_mode_(false),
        first_(first),
        last_(last),
        first2_(first2) {

  }
  void togle_two_sets_mode() {
    if (two_sets_)
      two_sets_mode_ = !two_sets_mode_;
  }
  void operator++() {
    if (first_ != last_) {
      ++first_;
      if (two_sets_)
        ++first2_;
    }
  }

  void render() {
    event_t event = *first_;
    track_t track = detector_->toPS(event);
    gp_->renderZYLine(
        track.z_dn(), -detector_->R(), track.z_up(), detector_->R());
    gp_->renderZYCircle(event.z(), event.y(), 4, true);

    if (two_sets_mode_) {
      gp_->set_color(glm::vec4(1, 0, 1, 1));
      event = *first2_;
      track = detector_->toPS(event);
      gp_->renderZYLine(
          track.z_dn(), -detector_->R(), track.z_up(), detector_->R());
      gp_->renderZYCircle(event.z(), event.y(), 4, true);
      gp_->set_color(glm::vec4(0, 0, 0, 1));
    }
  }
 private:
  GeometryPlot* gp_;
  detector_t* detector_;
  bool two_sets_;
  bool two_sets_mode_;
  const_iterator first_;
  const_iterator last_;
  const_iterator first2_;
};
