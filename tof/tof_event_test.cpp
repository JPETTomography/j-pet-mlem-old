#include<cmath>
#include<gtest/gtest.h>



#include"tof_event.h"
#include"tof_detector.h"

class tof_event_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    
  }


};


TEST_F(tof_event_test,set_test) {
  ToF_Event_2D<double> event(200.0,150.0,1.0);

  ASSERT_DOUBLE_EQ(200,event.z());
  ASSERT_DOUBLE_EQ(150.0,event.y());
  ASSERT_DOUBLE_EQ(1.0,event.tan());

}

TEST_F(tof_event_test,track_event_track_test) {
  ToF_Event_2D<double> event;
  ToF_Detector_2D<double> detector(350.0,500.0);

  ToF_Track_2D<double> track(350,-350,0.0);
  
  event=detector.fromPS(track);


  ASSERT_DOUBLE_EQ(0.0,event.z());
  ASSERT_DOUBLE_EQ(0.0,event.y());
  ASSERT_DOUBLE_EQ(1.0,event.tan());

  ToF_Track_2D<double> rec_track;

  rec_track=detector.toPS(event);

  ASSERT_DOUBLE_EQ(track.z_up(),rec_track.z_up());
  ASSERT_DOUBLE_EQ(track.z_dn(),rec_track.z_dn());
  ASSERT_DOUBLE_EQ(track.dl(),rec_track.dl());

}

TEST_F(tof_event_test,event_track_event_test) {
  ToF_Event_2D<double> event(100,150,tan(M_PI/3.0));
  ToF_Detector_2D<double> detector(350.0,500.0);

  ToF_Track_2D<double> track=detector.toPS(event);
  
  ToF_Event_2D<double> rec_event=detector.fromPS(track);

  ASSERT_DOUBLE_EQ(event.tan(),rec_event.tan());
  ASSERT_DOUBLE_EQ(event.y(),rec_event.y());
  ASSERT_DOUBLE_EQ(event.z(),rec_event.z());

}

int 
main(int argc,char *argv[]) {
::testing::InitGoogleTest(&argc,argv);

  return RUN_ALL_TESTS();
  
}
