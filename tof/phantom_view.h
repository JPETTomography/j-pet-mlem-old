#ifndef __PHANTOM_VIEW_H__
#define __PHANTOM_VIEW_H__

#include "geometry_plot.h"
#include "phantom.h"
class PhantomView {
public:
  PhantomView(GeometryPlot *gp, Phantom *phantom):gp_(gp), phantom_(phantom) {}

  void render() {
    Phantom::const_iterator it = phantom_->begin();
    for(;it! = phantom_->end();++it) {
      std::cerr << "next region " << std::endl;
      std::cerr << "x " << (*it)->x() << std::endl;
      gp_->renderZYEllipse((*it)->x(), (*it)->y(), (*it)->a(), (*it)->b(), (*it)->phi()
         );
      std::cerr << "end region " << std::endl;
    }
  }

private:
  GeometryPlot *gp_;
  Phantom *phantom_;
};

#endif
