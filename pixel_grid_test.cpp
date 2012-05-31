#include<cmath>
#include<gtest/gtest.h>



#include"pixel_grid.h"

class pixel_grid_test : public ::testing::Test {
protected:

  pixel_grid_test():
    ll(Point<>(-250,-200)),
    ur(Point<>(100,50)),
    nx(70),ny(125),grid(ll,ur,nx,ny) {};
			    
		    

  virtual void SetUp() {
    dx=5;
    dy=2;
  }

#if 1
  int nx;
  int ny;

  Point<> ll;
  Point<> ur;
  PixelGrid<> grid;

  double dx,dy;
#endif
  
};


TEST_F(pixel_grid_test,point_index_set_test) {
  double x=1.0;
  double y=2.8;
  Point<> p(x,y); 

  ASSERT_DOUBLE_EQ(x,p.x);
  ASSERT_DOUBLE_EQ(y,p.y);
  
  int ix=10;
  int iy=28;
  Index<> ind(ix,iy); 

  ASSERT_EQ(ix,ind.x);
  ASSERT_EQ(iy,ind.y);
  
}

TEST_F(pixel_grid_test,set_test) {

  ASSERT_EQ(nx,grid.nx());
  ASSERT_EQ(ny,grid.ny());
  
  ASSERT_DOUBLE_EQ(dx,grid.dx());
  ASSERT_DOUBLE_EQ(dy,grid.dy());
  
}


TEST_F(pixel_grid_test,index_test) {


  ASSERT_EQ(23*70+13,grid.index(13,23));

  
}
