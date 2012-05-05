#include<cmath>
#include<cstdio>
#include<cstdlib>

#include<vector>

#include<boost/program_options/options_description.hpp>
#include<boost/program_options/variables_map.hpp>
#include<boost/program_options/parsers.hpp>
namespace po = boost::program_options;

#include"event.h"
#include"detector.h"
#include"reconstruction.h"

main(int argc, char *argv[]) {

  po::options_description desc("Allowed options");
  desc.add_options()
    ("n-pixels,n", po::value<int>(), "set number of  rind scanners")
    ("n-detectors,d", po::value<int>(), "set number of pixels in one dimension");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm); 

  if(vm.count("n-detectors")==0 || vm.count("n-pixels")==0) {
    std::cerr<<" number of detectors or number of pixels not set"<<std::endl;
    exit(-1);
  } 
  std::cout<<vm["n-detectors"].as<int>()<<" ";
  std::cout<<vm["n-pixels"].as<int>()<<std::endl;


  Ring2DDetector<int,double> detector(vm["n-detectors"].as<int>(),
				      vm["n-pixels"].as<int>());


  ProbabilityMatrix<int,double> probability_matrix(detector);

  const int max_lors=detector.n_detectors()*detector.n_detectors();

  double *p = new double[max_lors];

  
  
  for(int ipix=0;ipix<probability_matrix.octant_size();++ipix) {
    Pixel<int,double> pix=probability_matrix.octant(ipix);

    for(int i=0;i<max_lors;++i)
      p[i]=0.0;

    fill_pixel(pix,detector.n_detectors(),p,1000000);

    std::vector<Lor<int,double> > *lrs=new std::vector<Lor<int,double> >;
    double sum=0.0;
    for(int i=0;i<max_lors;++i) {
      if(p[i]>0.0) {
	int f=i/detector.n_detectors();
	int s=i%detector.n_detectors();
	lrs->push_back(Lor<int,double>(f,s,p[i]));
	sum+=p[i];
      }
    }

    Row<Pixel<int,double>, double> pixel_row(pix,*lrs);
    delete lrs;

    std::cout<<pixel_row;


  }


}


