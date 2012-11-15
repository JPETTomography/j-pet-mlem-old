// PET Phantom
// Authors:
//   Piotr Bialas    <piotr.bialas@uj.edu.pl>
//   Jakub Kowal     <jakub.kowal@uj.edu.pl>
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//
// Generates phantom measurements using Monte Carlo.

#include <random>
#include <iostream>


#include <cmdline.h>
#include "cmdline_types.h"

#include"phantom.h"
#include "detector_ring.h"
#include "model.h"


#if _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {

try {
  cmdline::parser cl;
  cl.footer("pixels_file measurements_file");

#if _OPENMP
  cl.add<size_t> ("n-threads",   't', "number of OpenMP threads",false);
#endif
  cl.add<size_t>     ("n-pixels",    'n', "number of pixels in one dimension", false, 256);
  cl.add<size_t>     ("n-detectors", 'd', "number of ring detectors",          false, 64);
  cl.add<size_t>     ("n-emissions", 'e', "emissions",               false, 1);
  cl.add<double>     ("radious",     'r', "inner detector ring radious",       false);
  cl.add<double>     ("s-pixel",     'p', "pixel size",                        false);
  cl.add<double>     ("w-detector",  'w', "detector width",                    false);
  cl.add<double>     ("h-detector",  'h', "detector height",                   false);
  cl.add<std::string>("model",       'm', "acceptance model",                  false,
                      "scintilator", cmdline::oneof<std::string>("always", "scintilator"));
  cl.add<double>     ("acceptance",  'a', "acceptance probability factor",     false, 10.);
  cl.add             ("stats",         0, "show stats");
  cl.add<std::string>("output",      'o', "output binary triangular sparse system matrix", false);
  cl.add<std::mt19937::result_type>
    ("seed",        's', "random number generator seed",      false);
  cl.add<std::string>("phantom",'f',"phantom description file",false,"");
  cl.parse_check(argc, argv);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<size_t>("n-threads"));
  }
#endif
  

  auto n_pixels    = cl.get<size_t>("n-pixels");
  auto n_detectors = cl.get<size_t>("n-detectors");
  auto n_emissions = cl.get<size_t>("n-emissions");
  auto radious     = cl.get<double>("radious");
  auto s_pixel     = cl.get<double>("s-pixel");
  auto w_detector  = cl.get<double>("w-detector");
  auto h_detector  = cl.get<double>("h-detector");
  auto acceptance  = cl.get<double>("acceptance");

  // automatic radious
  if (!cl.exist("s-pixel")) {
    if (!cl.exist("radious")) {
      s_pixel = 2. / n_pixels; // exact result
    } else {
      s_pixel = M_SQRT2 * radious / n_pixels;
    }
    std::cerr << "--s-pixel=" << s_pixel << std::endl;
  }

  // automatic detector size
  if (!cl.exist("w-detector")) {
    w_detector = 2 * M_PI * .9 * radious / n_detectors;
    std::cerr << "--w-detector=" << w_detector << std::endl;
  }
  if (!cl.exist("h-detector")) {
    h_detector = w_detector;
    std::cerr << "--h-detector=" << h_detector << std::endl;
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  if (cl.exist("seed")) {
    gen.seed(cl.get<std::mt19937::result_type>("seed"));
  }


  size_t tubes[n_detectors][n_detectors];
  for(int i=0;i<n_detectors;i++)
    for(int j=0;j<n_detectors;j++)
      tubes[i][j]=0;

  detector_ring<double> dr(n_detectors, n_pixels, s_pixel, radious, w_detector, h_detector);

  int n_emitted=0;

  std::uniform_real_distribution<> one_dis(0., 1.);
  std::uniform_real_distribution<> fov_dis(-dr.fov_radius(), dr.fov_radius());
  std::uniform_real_distribution<> phi_dis(0., M_PI);
 
  double fov_r2=dr.fov_radius()*dr.fov_radius();

  Phantom  phantom;

  if(cl.exist("phantom")) {
    std::string name=cl.get<std::string>("phantom");
    FILE *fin=fopen(name.c_str(),"r");
    if(fin==NULL) {
      fprintf(stderr,"cannot open file `%s' for reading\n",name.c_str());
      exit(-1);
    }
    phantom.load_from_file(fin);
  }
    ;

#if 0
  phantom.addRegion(0,0,0.2,0.25,0.0,0.1);
  phantom.addRegion(0,.1,0.025,0.1,0.0,0.7);
#endif

  scintilator_accept<> model(acceptance);

  while(n_emitted<n_emissions) {

    double x=fov_dis(gen);
    double y=fov_dis(gen);
    if(x*x+y*y < fov_r2) {

      if(phantom.emit(x,y,one_dis(gen))) {
        detector_ring<double>::lor_type lor; 
        std::cerr<<n_emitted<<" "<<x<<" "<<y<<"\n";
        double angle=phi_dis(gen);
        auto hits=dr.emit_event(gen,model,x,y,angle,lor);
        if(hits==2) {
          if(lor.first>lor.second)
            std::swap(lor.first,lor.second);
          tubes[lor.first][lor.second]++;
          // std::cerr<<"hits\n";
        }
        n_emitted++;
      }
    }
  }

#if 0
  for(int i=0;i<n_detectors;i++)
    for(int j=i+1;j<n_detectors;j++) {
            if(tubes[i][j]>0) 
        std::cout<<i<<" "<<j<<"  "<<tubes[i][j]<<"\n";
    }
#endif


  return 0;

 } catch(std::string &ex) {
  std::cerr << "error: " << ex << std::endl;
} catch(const char *ex) {
  std::cerr << "error: " << ex << std::endl;
}

}
