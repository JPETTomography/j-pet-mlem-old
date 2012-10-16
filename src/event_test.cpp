#include<cmath>
#include<gtest/gtest.h>


#include"event.h"
#include"detector.h"

class eventTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    
  }


};


TEST_F(eventTest,set_test) {
  event<double> event(1.0,0.5,2.0);

  ASSERT_DOUBLE_EQ(1.0,event.x());
  ASSERT_DOUBLE_EQ(0.5,event.y());
  ASSERT_DOUBLE_EQ(2.0,event.phi());

}

TEST_F(eventTest,detector_center_test) {
  event<double> event(0.0,0.0,2.0);

  std::pair<double,double> time=tof(event,1.0);
  
  ASSERT_DOUBLE_EQ( 1.0,time.first);
  ASSERT_DOUBLE_EQ(-1.0,time.second);
  


}

TEST_F(eventTest,detector_on_x_axis_test) {
  event<double> event1(0.5,0.0,0.0);

  std::pair<double,double> time=tof(event1,1.0);
  
  ASSERT_DOUBLE_EQ( 0.5,time.first);
  ASSERT_DOUBLE_EQ(-1.5,time.second);
  
  event<double> event2(0.5,0.0,M_PI);

  time=tof(event2,1.0);
  
  ASSERT_DOUBLE_EQ( 1.5,time.first);
  ASSERT_DOUBLE_EQ(-0.5,time.second);

}

TEST_F(eventTest,detector_on_y_axis_test) {
  event<double> event1(0.0,0.5,M_PI/2.0);

  std::pair<double,double> time=tof(event1,1.0);
  
  ASSERT_DOUBLE_EQ( 0.5,time.first);
  ASSERT_DOUBLE_EQ(-1.5,time.second);
  
  event<double> event2(0.0,0.5,-M_PI/2.0);

  time=tof(event2,1.0);
  
  ASSERT_DOUBLE_EQ( 1.5,time.first);
  ASSERT_DOUBLE_EQ(-0.5,time.second);

}


TEST_F(eventTest,detector_on_sym_test) {
  double x=0.5;
  double y=0.5;

  event<double> event1(x,y,M_PI/4.0);

  std::pair<double,double> time=tof(event1,1.0);
  
  double dist=sqrt(x*x+y*y);
  
  ASSERT_DOUBLE_EQ( 1.0 - dist,time.first);
  ASSERT_DOUBLE_EQ(-1.0 - dist,time.second);
  

  event<double> event2(x,y,5.0*M_PI/4.0);

  time=tof(event2,1.0);
  
  ASSERT_DOUBLE_EQ( 1.0 + dist,time.first);
  ASSERT_DOUBLE_EQ(-1.0 + dist,time.second);
  
}

TEST_F(eventTest,detector_radius_test) {
  double x= 0.8;
  double y=-0.6;

  event<double> event(x,y,2.1);

  std::pair<double,double> time=tof(event,1.0);
  
  double s,c;
  sincos(event.phi(),&s,&c);
  
  double x_hit=x+time.first*c;
  double y_hit=y+time.first*s;

  ASSERT_DOUBLE_EQ(1.0,x_hit*x_hit+y_hit*y_hit);

  x_hit=x+time.second*c;
  y_hit=y+time.second*s;

  ASSERT_DOUBLE_EQ(1.0,x_hit*x_hit+y_hit*y_hit);


}


TEST_F(eventTest,detector_lor_center_test) {
  event<double> event1(0.0,0.0,0.001);

  std::pair<double,double> time=tof(event1,1.0);
  std::pair<short,short> lors=lor(time,event1,1.0,128);
  
  ASSERT_EQ(0,lors.first);
  ASSERT_EQ(64,lors.second);

  event<double> event2(0.0,0.0,M_PI/2.0+0.001);

  time=tof(event2,1.0);
  lors=lor(time,event2,1.0,128);
  
  ASSERT_EQ(32,lors.first);
  ASSERT_EQ(96,lors.second);

  event<double> event3(0.0,0.0,3.0*M_PI/2.0+0.001);

  time=tof(event3,1.0);
  lors=lor(time,event3,1.0,128);
  
  ASSERT_EQ(32,lors.first);
  ASSERT_EQ(96,lors.second);


}

int 
main(int argc,char *argv[]) {
::testing::InitGoogleTest(&argc,argv);

  return RUN_ALL_TESTS();
  
}
