#include<cmath>
#include<gtest/gtest.h>



#include"tof_event.h"
#include"tof_detector.h"
#include"phantom.h"


class elliptical_region_test : public ::testing::Test {
protected:
 
  virtual void SetUp() {    
    disk=new EllipticalRegion(1.0,1.0,2.0,2.0,0.0,0.5);
    region= new EllipticalRegion(0,1,1,0.5,M_PI/3.0,0.75); 
   
   

   
  }
 

  EllipticalRegion *disk;
  EllipticalRegion *region;
  
  
};


TEST_F(elliptical_region_test, getter_test) {
  ASSERT_DOUBLE_EQ(disk->activity(),0.5);
  ASSERT_DOUBLE_EQ(region->activity(),0.75);
}


TEST_F(elliptical_region_test, in_test) {
  ASSERT_EQ(true,disk->in(1,1));
  ASSERT_EQ(true,disk->in(1.563,-0.8545));
  ASSERT_EQ(false,disk->in(-0.677,-2.5));

  ASSERT_EQ(true,region->in(-0.328,0.26));
  ASSERT_EQ(true,region->in(0.4371,1.792));
  ASSERT_EQ(false,region->in(1,1));
  
}



class phantom_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    //    phantom = new Phantom(-3,-3,3,3);
    phantom = new Phantom;
    phantom->addRegion(1.0,1.0,2.0,2.0,0.0,0.5);
    phantom->addRegion(0,1,1,0.5,M_PI/3.0,0.75); 
  }
  
  Phantom *phantom;
  
};
  
TEST_F(phantom_test,activity_test) {

  ASSERT_DOUBLE_EQ(0.5 , phantom->activity(1,1));
  ASSERT_DOUBLE_EQ(0.5 , phantom->activity(1.563,-0.8545));
  ASSERT_DOUBLE_EQ(0.0 , phantom->activity(-0.677,-2.5));

  ASSERT_DOUBLE_EQ(0.75 , phantom->activity(-0.328,0.26));
  ASSERT_DOUBLE_EQ(0.75 , phantom->activity(0.4371,1.792));
  
}


TEST_F(phantom_test,emit_test) {

  ASSERT_EQ(false , phantom->emit(1,1,.75));
  ASSERT_EQ(true , phantom->emit(1,1,.45));
  ASSERT_EQ(false , phantom->emit(1.563,-0.8545,.51));
  ASSERT_EQ(true , phantom->emit(1.563,-0.8545,.10));
  ASSERT_EQ(false, phantom->emit(-0.677,-2.5,.1));
  ASSERT_EQ(false, phantom->emit(-0.677,-2.5,.25));
  ASSERT_EQ(false, phantom->emit(-0.677,-2.5,0.001));  

  ASSERT_EQ(false , phantom->emit(-0.328,0.26,.76));
  ASSERT_EQ(true , phantom->emit(-0.328,0.26,.74));
  ASSERT_EQ(false , phantom->emit(0.4371,1.792,.77));
  ASSERT_EQ(true , phantom->emit(0.4371,1.792,.73));
  
}
