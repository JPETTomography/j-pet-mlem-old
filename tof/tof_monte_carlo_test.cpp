#include <cmath>
#include <gtest/gtest.h>
#include "tof_event.h"
#include "tof_detector.h"
#include "phantom.h"
#include "tof_monte_carlo.h"

class tof_monte_carlo_test : public ::testing::Test {
protected:
  tof_monte_carlo_test():detector(350, 500) {};
  virtual void SetUp() {
    detector.set_sigma(11, 32);
    mc = new   ToF_Monte_Carlo<double, detector_t>(detector);
    mc->gen_seeds(544445);

  }
  typedef ToF_Detector_2D<double> detector_t;
  typedef ToF_Event_2D<double> event_t;
  detector_t    detector;
  ToF_Monte_Carlo<double, detector_t> *mc;
};
TEST_F(tof_monte_carlo_test, add_noise_to_track_test) {
  ToF_Track_2D<double> track_in(-200, 23, 117);

  ToF_Track_2D<double> track_out = mc->add_noise(track_in);

  ASSERT_NE(track_in.z_up(), track_out.z_up());
  ASSERT_NE(track_in.z_dn(), track_out.z_dn());
  ASSERT_NE(track_in.dl(), track_out.dl());
}

TEST_F(tof_monte_carlo_test, add_noise_to_event_test) {
  ToF_Event_2D<double> event_in(-10, 23, -.75);

  ToF_Event_2D<double> event_out = mc->add_noise(event_in);

  ASSERT_NE(event_in.z(), event_out.z());
  ASSERT_NE(event_in.y(), event_out.y());
  ASSERT_NE(event_in.tan(), event_out.tan());
}

