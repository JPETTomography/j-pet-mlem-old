#ifndef __TOF_EVENT__
#define __TOF_EVENT__

template<typename F> class  ToF_Track_2D {

public:
  typedef F float_t;

  ToF_Track_2D(){};

  ToF_Track_2D(F z_up,F z_dn, F dl): z_up_(z_up), z_dn_(z_dn), dl_(dl){};

  F z_up() const {return z_up_;}
  F z_dn() const {return z_dn_;}
  F dl()   const {return dl_;}
  
private:
  F z_up_;
  F z_dn_;
  F dl_;
};

template<typename F> class  ToF_Event_2D {


public:
  typedef F float_t;


  ToF_Event_2D(){};
  ToF_Event_2D(F z, F y, F tan):z_(z),y_(y),tan_(tan){}
    
  
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
