#ifndef __TOPET_SIMULATOR_H__
#define __TOPET_SIMULATOR_H__


#include"tof_monte_carlo.h"
#include"pixel_grid.h"


template<typename F> class TOPETSimulator {

public:

  typedef ToF_Event_2D<F> event_t;
  typedef ToF_Track_2D<F> track_t;
  

  void init() {
  }


  void simulate() {
  }

  

private:
  
  ToF_Monte_Carlo *mc_;
  
  PixelGrid<F> emitted_density_;
  PixelGrid<F> detected_density_;
  PixelGrid<F> tof_density_;
  
  std::vector<event_t> emitted_events_;
  std::vector<event_t> detected_events_;
  std::vector<event_t> tof_events_;

};



#endif
