#ifndef __TOPET_SIMULATOR_H__
#define __TOPET_SIMULATOR_H__
#include "tof_monte_carlo.h"
#include "pixel_grid.h"
template<typename F> class TOPETSimulator {

public:

  typedef ToF_Event_2D<F> event_t;
  typedef ToF_Track_2D<F> track_t;
  typedef Point<F> point_t;
  typedef ToF_Detector_2D<F> detector_t;
  typedef typename std::vector<event_t>::const_iterator const_iterator;

  void init() {

    int width = 100;
    int height = 100;
    emitted_density_ = new PixelGrid<F>(point_t(-250, -250), point_t(250, 250), width, height);
    detected_density_ = new PixelGrid<F>(point_t(-250, -250), point_t(250, 250), width, height);
    tof_density_ = new PixelGrid<F>(point_t(-250, -250), point_t(250, 250), width, height);
    acceptance_density_ = new PixelGrid<F>(point_t(-250, -250), point_t(250, 250), width, height);

    detector_ = new detector_t(350, 500);
    detector_->set_sigma(11, 32);

    mc_ = new ToF_Monte_Carlo<F, detector_t>(*detector_);
    mc_->gen_seeds(5565665);
    phantom_ = new Phantom;
    phantom_->addRegion(0.0, 0.0, 200.0, 100.0, 0.0, 0.25);
    phantom_->addRegion(-100.0, 0.0, 50, 70.0, M_PI/3.0, 0.50);
    phantom_->addRegion(150, 20, 10, 17, M_PI/4.0, 0.80);
    fov_ll_x_ = -250.0;
    fov_ll_y_ = -250.0;

    fov_ur_x_ = 250.0;
    fov_ur_y_ = 250.0;
    for(int ix = 0;ix<acceptance_density_->nx();ix++)
      for(int iy = 0;iy<acceptance_density_->ny();iy++) {
  point_t c = acceptance_density_->center(ix, iy);
    F acc = detector_->acceptance(c.x, c.y);
  acceptance_density_->add(ix, iy, acc);
      }

    F max = acceptance_density_->max();
    acceptance_density_->divide_by(max);

  }

   detector_t *detector()  {return detector_;}
   Phantom *phantom()  {return phantom_;}

  void simulate_from_phantom(int n) {
    std::cerr << "simulating " << n << " events" << std::endl;
    emitted_events_.clear();
    emitted_events_.resize(n);
    int n_emitted_events = mc_->fill_with_events_from_phantom(emitted_events_.begin(), phantom_, fov_ll_x_, fov_ll_y_, fov_ur_x_, fov_ur_y_, n
              );

    std::cerr << "emited " << n_emitted_events << " events" << std::endl;

    std::cerr << "inserting " << n << " events" << std::endl;
    emitted_density_->insert(emitted_events_.begin(), n_emitted_events);
    std::cerr << "inserted " << n << " events" << std::endl;

    F max = emitted_density_->max();
    emitted_density_->divide_by(max);

    detected_events_.clear();
    detected_events_.resize(n_emitted_events);

    std::cerr << "adding  noise " << std::endl;
    typename std::vector<event_t>::iterator it = detected_events_.begin();
    int n_detected_events = 0;
    for(int i = 0;i<n_emitted_events;i++) {
      if(detector_->detected(emitted_events_[i])) {
  (*it) = emitted_events_[i];
  ++it;
  n_detected_events++;
      }
    }
    detected_density_->insert(detected_events_.begin(), n_detected_events);
    F detected_max = detected_density_->max();
    detected_density_->divide_by(max);

    tof_events_.clear();
    tof_events_.resize(n_detected_events);

    std::cerr << "adding  noise " << std::endl;
    mc_->add_noise(
       detected_events_.begin(), detected_events_.begin()+n_detected_events, tof_events_.begin()
       );
    std::cerr << "added   noise to " << n_detected_events << std::endl;

    std::cerr << "inserting " << n_detected_events << " events in tof" << std::endl;
    tof_density_->insert(tof_events_.begin(), n_detected_events);
    std::cerr << "inserted " << n_detected_events << " events in tof" << std::endl;

    F tof_max = tof_density_->max();
    tof_density_->divide_by(tof_max);
  }

  void simulate_from_single_point(int n) {
    std::cerr << "simulating " << n << " events" << std::endl;
    emitted_events_.clear();
    emitted_events_.resize(n);

    std::cerr << "emiting " << n << " events" << std::endl;
    int n_emitted_events = mc_->fill_with_events_from_single_point(emitted_events_.begin(), 0.0f, 0.0f, n);
    std::cerr << "emited " << n_emitted_events << " events" << std::endl;

    std::cerr << "inserting " << n << " events" << std::endl;
    emitted_density_->insert(emitted_events_.begin(), n_emitted_events);
    std::cerr << "inserted " << n << " events" << std::endl;

    F max = emitted_density_->max();
    emitted_density_->divide_by(max);

    tof_events_.clear();
    tof_events_.resize(n_emitted_events);

    std::cerr << "adding  noise " << std::endl;
    int n_detected_events = mc_->add_noise_to_detected(
         emitted_events_.begin(), emitted_events_.begin()+n_emitted_events, tof_events_.begin()
         );
    std::cerr << "added   noise to " << n_detected_events << std::endl;

    std::cerr << "inserting " << n_detected_events << " events in tof" << std::endl;
    tof_density_->insert(tof_events_.begin(), n_detected_events);
    std::cerr << "inserted " << n_detected_events << " events in tof" << std::endl;

    F tof_max = tof_density_->max();
    tof_density_->divide_by(tof_max);
  }
  const  F *emitted_density() {return emitted_density_->pixels();}
  const  F *detected_density() {return detected_density_->pixels();}
  const  F *tof_density() {return tof_density_->pixels();}
  const  F *acceptance_density() {return acceptance_density_->pixels();}
  const_iterator emitted_begin() {return emitted_events_.begin();}
  const_iterator emitted_end() {return emitted_events_.end();}

  const_iterator detected_begin() {return detected_events_.begin();}
  const_iterator detected_end() {return detected_events_.end();}

  const_iterator tof_begin() {return tof_events_.begin();}
  const_iterator tof_end() {return tof_events_.end();}

  ~TOPETSimulator() {
    delete detector_;
    delete phantom_;
  }

private:
  F fov_ll_x_;
  F fov_ll_y_;
  F fov_ur_x_;
  F fov_ur_y_;

  detector_t *detector_;
  Phantom *phantom_;

  ToF_Monte_Carlo<F, detector_t> *mc_;

  PixelGrid<F> *emitted_density_;
  PixelGrid<F> *detected_density_;
  PixelGrid<F> *tof_density_;
  PixelGrid<F> *acceptance_density_;

  std::vector<event_t> emitted_events_;
  std::vector<event_t> detected_events_;
  std::vector<event_t> tof_events_;
};
#endif