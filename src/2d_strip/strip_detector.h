#include"event.h"
/**
 * Class  responsible for the Strip detector together with
 * the pixel grid inside. 
 */
template<typename F> class StripDetector {
public:
  StripDetector(F radius,
                F scintilator_length,
                int n_y_pixels,
                int n_z_pixels,
                F pixel_height, //y direction
                F pixel_width,  //z direction
                F sigma_z,
                F sigma_dl,
                F grid_center_y = 0.0,
                F grid_center_z = 0.0):
    radius_(radius),
    scintilator_length_(scintilator_length),
    n_y_pixels_(n_y_pixels),
    n_z_pixels_(n_z_pixels), 
    pixel_width_(pixel_width),
    pixel_height_(pixel_height),
    sigma_z_(sigma_z),
    sigma_dl_(sigma_dl),
    grid_center_y_(grid_center_y),
    grid_center_z_(grid_center_z) {    
  }


  F radius() const {return radius_;}
  F scinatilator_length() const {return scintilator_length_;}
  F half_scintilator_length() const {return F(0.5)*scintilator_length_;}
  F s_z()  const {return sigma_z_;}
  F s_dl() const {return sigma_dl_;}

  event<F> to_projection_space_tan(const ImageSpaceEventTan<F>& ev) {    
    F z_u = ev.z + (radius()- ev.y)*ev.tan;
    F z_d = ev.z - (radius()+ ev.y)*ev.tan;
    F dl  = -F(2.0)* ev.y*sqrt(ev.tan*ev.tan+F(1.0));
    return event<F>(z_u, z_d, dl);
  }

  event<F> to_projection_space_angle(const ImageSpaceEventAngle<F>& img_event) {
    return to_angle(to_projection_space_angle(img_event));
  }

  ImageSpaceEventTan<F> from_projection_space_tan(const event<F> &ev) {
    F t = event_tan(ev.z_u, ev.z_d, radius() );
    F y = event_y(ev.dl, t);
    F z = event_z(ev.z_u, ev.z_d, y, t);
    return ImageSpaceEventTan<F>(y, z, t);
  }

  ImageSpaceEventAngle<F> from_projection_space_angle(const event<F> &ev) {
    return to_angle(from_projection_space_tan(ev));
  }

private:
  const F radius_;
  const F scintilator_length_;
  const int n_y_pixels_;
  const int n_z_pixels_;
  const F pixel_width_;
  const F pixel_height_;
  const F sigma_z_;
  const F sigma_dl_;
  const F grid_center_y_;
  const F grid_center_z_;
  
};
