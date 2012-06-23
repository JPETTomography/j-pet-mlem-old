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

  int n_emissions;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("n-pixels,n", po::value<int>(), "set number of  rind scanners")
    ("n-detectors,d", po::value<int>(), "set number of pixels in one dimension")    ("n-emissions,e", po::value<int>(&n_emissions)->default_value(1000000), "emissions per pixel");

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

    fill_pixel(pix,detector.n_detectors(),p,n_emissions);


    Row<Pixel<int,double>, double> *pixel_row=row_from_array(pix,p,detector.n_detectors(),max_lors);

    probability_matrix.push_back_row(pixel_row);
    std::cout<<*pixel_row;


  }

  
  FILE *fout=fopen("prob.bin","w");
  probability_matrix.fwrite(fout);
  fclose(fout);

}


