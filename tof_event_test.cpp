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

TEST_F(tof_event_test,formPS_test) {
  ToF_Event_2D<double> event;
  ToF_Detector_2D<double> detector(350.0,500.0);

  event=detector.fromPS(350.0,-350.0,0.0);


  ASSERT_DOUBLE_EQ(0.0,event.z());
  ASSERT_DOUBLE_EQ(0.0,event.y());
  ASSERT_DOUBLE_EQ(1.0,event.tan());

}

int 
main(int argc,char *argv[]) {
::testing::InitGoogleTest(&argc,argv);

  return RUN_ALL_TESTS();
  
}
