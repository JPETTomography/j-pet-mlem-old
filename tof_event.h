#ifndef __TOF_EVENT__

template<typename F> class  ToF_Event_2D {


public:
  typedef F float_t;

  static void fromPS(ToF_Event_2D &event,F z_up, F z_dn, F dl,F R) {
    event.tan_=(z_up-z_dn)/((F)2.0*R);
    event.z_=(F)(0.5)*(z_up+z_dn);
    event.y_=-dl/sqrt(event.tan_*event.tan_+1);
  }

  ToF_Event_2D(){};
  ToF_Event_2D(F z, F y, F t):z_(z),y_(y),tan_(t){}
    
  
  F z() const {return z_;}
  F y() const {return y_;}
  F tan() const {return tan_;}
  

private:
  
  F z_;
  F y_;
  F tan_;

};


struct row {
  int r;
  int begin;
  int end;
};

template<typename F> class Event_Pixels {

private:
  int pixel_;
  std::vector<row > rows_;
  std::vector<int>  pixels_;
};



template<typename F> class pixel_grid {
public:
  typedef F float_t;
  
private:
  int n_z_;
  int n_y_;
  
  F ll_z_;
  F ll_y_;
  F ur_z_;
  F ur_y_;
  
  
};


#endif
