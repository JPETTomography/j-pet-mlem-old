#pragma once
#pragma once

#include <vector>
#include <list>

#define LOCATION(x,y,size) (x*size) + y
#define CRICLE_REGION(x,y,r) ((x-(r/2))*(x-(r/2))) + ((y-(r/2))*(y-(r/2))) < (((r+1)/2)*((r+1)/2)) ? true : false
#define FALSE_CRICLE_REGION(x,y,r) ((x*x) + (y*y)) > (r*r) ? true : false

typedef std::pair<size_t, size_t> pixel_location;
typedef uint32_t file_int;
typedef uint16_t file_half;

 
 #define fourcc(a, b, c, d) (((d)<<24) | ((c)<<16) | ((b)<<8) | (a))
 
// binary serialization                                  // n_pixels  n_detectors  triagular
static const file_int magic_1 = fourcc('P','E','T','t'); //                           X
static const file_int magic_2 = fourcc('P','E','T','s'); //     X                     X
static const file_int magic_t = fourcc('P','E','T','p'); //     X          X          X
static const file_int magic_f = fourcc('P','E','T','P'); //     X          X


using namespace std;

struct hits_per_pixel{

	pixel_location location;
	float probability;

};


class Pet_reconstruction{

private:
	int n_tubes;
	int n_pixels;
	int n_iter;
	int emissions;
	int n_lors;
	static vector< vector<hits_per_pixel> > system_matrix;
	static vector<float> n;
	float *scale;


public:
	Pet_reconstruction(int n_iter,string file):
	n_iter(n_iter)
	{	
		
		ibstream in(file);
		
		typedef uint32_t file_int;
		typedef uint16_t file_half;
		
	  
		//ibstream in("pet10k_full.bin");
		
		file_int in_magic;
		
		in >> in_magic;
		
		in >> n_pixels;
		in >> emissions;
		in >> n_tubes;	
		
		printf("EM: %d\n",emissions);
		
		scale = (float*)malloc(n_pixels*n_pixels*sizeof(float));
		memset(scale,0.1f,sizeof(scale));
		
		printf("%f %f %f\n",scale[0],scale[1],scale[2]);

		std::vector<hits_per_pixel> pixels;
		
		
		
		for(;;){
	  
			if (in.eof()) break;
			
			file_half a, b;
			in >> a >> b;

			file_int count;

			in >> count;
				
				int sum = 0;
				
				//cout << "blabla" << endl;
				for(int i = 0;i < count; ++i){
					
					file_half x, y;
					file_int hits;
		
					in >> x >> y >> hits;
					pixel_location pixel(x,y);
					hits_per_pixel data;
					data.location = pixel;
					data.probability = static_cast<float>(hits/static_cast<float>(emissions));
					
					scale[(x*n_pixels) + y] += data.probability;  
				
					pixels.push_back(data);
			  
					sum+=hits;
	
				}
				
			n.push_back(sum);
			
			system_matrix.push_back(pixels);
			pixels.clear();
	
		}
		
		printf("Pixels: %d\nEmissions: %d\nTubes: %d\nLors: %d\n",n_pixels,emissions,n_tubes,(int)system_matrix.size());
	
	}

	~Pet_reconstruction(){}

	std::vector<float> emt(){
		
		cout << n_tubes << " " << n_pixels << endl; 	
		
		float y[n_pixels * n_pixels];
		float u[n_lors];
		std::vector<float> rho(n_pixels * n_pixels,0.5);
		
		std::vector< vector<hits_per_pixel> >::iterator it;
		vector<hits_per_pixel>::iterator lt;
		
		for(int i = 0; i < n_iter;++i){
		
			for(int p = 0; p < n_pixels * n_pixels;++p){

				y[p] = 0.0;
			
			}

				int t = 0;
				
				for(it=system_matrix.begin(); it != system_matrix.end(); ++it){
					
					u[t] = 0.f;
					
					for(lt = (*it).begin(); lt != (*it).end();++lt){
					
						u[t] += (rho[LOCATION((*lt).location.first,(*lt).location.second,n_pixels)] * (*lt).probability) ;

					}

					t++;
									
				}
				
				t = 0;
	
				for(it=system_matrix.begin(); it != system_matrix.end(); ++it){
					
					
					float phi = n[t]/u[t];
					
					for(lt = (*it).begin(); lt != (*it).end();++lt){
					
						if(CRICLE_REGION((*lt).location.first,(*lt).location.second,n_pixels)){
						
							y[LOCATION((*lt).location.first,(*lt).location.second,n_pixels)] += phi * (*lt).probability;
						
						}
					}
					t++;
									
				}
			int valCond = 0;
			for(int width = 0; width < n_pixels; ++width)
			{
				for(int height = 0; height < n_pixels; ++height)
				{
			
					if(CRICLE_REGION(height,width,n_pixels))
					{
						rho[LOCATION(width,height,n_pixels)] = (rho[LOCATION(width,height,n_pixels)]/scale[LOCATION(width,height,n_pixels)]) * y[LOCATION(width,height,n_pixels)];
						 //scale[LOCATION(width,height,n_pixels)]
						//printf("%d %d %d rho: %f\n",width,height,width * n_pixels + height,rho[LOCATION(width,height,n_pixels)]);
						valCond++;
					}
				}
			}
			printf("Non Zero: %d\n",valCond);
		}
		
		return rho;
		
	}
	
	int get_n_pixels(){return n_pixels;}
	
	
};


