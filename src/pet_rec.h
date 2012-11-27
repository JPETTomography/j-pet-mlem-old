#pragma once

#include<iostream>
#include <vector>
#include <list>

#define LOCATION(x,y,size)  (x*size) +y

typedef std::pair<size_t, size_t> pixel_location;
typedef std::pair<size_t, size_t> lor;
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

struct mean_per_lor{

  lor Lor;
  int n;

};


class Pet_reconstruction{

public:
  Pet_reconstruction(int n_iter,string matrix,string mean):
  n_iter(n_iter)
  {

    ibstream in(matrix);

    typedef uint32_t file_int;
    typedef uint16_t file_half;

    file_int in_magic;

    in >> in_magic;
    in >> n_pixels;
    in >> emissions;
    in >> n_tubes;

    scale = (float*)malloc(n_pixels*n_pixels*sizeof(float));
    for(int p = 0; p < n_pixels * n_pixels; ++p){

      scale[p] = 0.f;

    }

    std::ifstream mean_file(mean);
    std::vector<mean_per_lor> lor_mean;


    for(;;)
    {

      int x,y,value;

      if (mean_file.eof()) break;

      mean_file >> y >> x >> value;

      mean_per_lor temp_obj;

      lor tmp(x,y);

      temp_obj.Lor = tmp ;
      temp_obj.n = value;

      lor_mean.push_back(temp_obj);

    }


    std::vector<hits_per_pixel> pixels;

    int index = 0;

    for(;;)
    {

      if (in.eof()) break;

      file_half a, b;
      in >> a >> b;
      lor _lor(a,b);

      n.push_back(get_mean_per_lor(a,b,lor_mean));

      //if(get_mean_per_lor(a,b,lor_mean) != 0){printf("get_mean: %d\n",get_mean_per_lor(a,b,lor_mean));}

      file_int count;

      in >> count;

        //cout << "blabla" << endl;
        for(int i = 0;i < count; ++i)
        {

          file_half x, y;
          file_int hits;

          in >> x >> y >> hits;
          pixel_location pixel(x,y);
          hits_per_pixel data;
          data.location = pixel;
          data.probability = static_cast<float>(hits/static_cast<float>(emissions));

          scale[LOCATION(x,y,n_pixels)] += data.probability;

          pixels.push_back(data);

        }

      system_matrix.push_back(pixels);
      pixels.clear();
      index++;
    }

    printf("Pixels: %d\nEmissions: %d\nTubes: %d\nLors: %d\n",n_pixels,emissions,n_tubes,(int)system_matrix.size());
    printf("%d\n",system_matrix.size());

  }

  ~Pet_reconstruction(){free(scale);}

  std::vector<float> emt(){

    float y[n_pixels * n_pixels];
    std::vector<float> u(system_matrix.size(),0.f);
    std::vector<float> rho(n_pixels * n_pixels,1.0f);

    for(int i = 0; i < n_iter;++i)
    {
      if(i%10){
        printf("%d\n",i);
      }

      for(int p = 0; p < n_pixels * n_pixels;++p)
      {

        y[p] = 0.f;

      }

      vector< vector<hits_per_pixel> >::iterator it_vector;
      vector<hits_per_pixel>::iterator it_list;

       int t = 0;

                for(it_vector = system_matrix.begin(); it_vector != system_matrix.end(); it_vector++)
          {
            u[t] = 0.0f;

                        for(it_list = (*it_vector).begin(); it_list != (*it_vector).end(); it_list++)
                        {

              u[t] += rho[LOCATION((*it_list).location.first,(*it_list).location.second,n_pixels)] * (*it_list).probability ;

                        }
            t++;

                    }

        t = 0;

        for(it_vector=system_matrix.begin(); it_vector != system_matrix.end(); ++it_vector)
        {

          if(n[t]> 0){

            float phi = n[t]/u[t];// n[t]/u[t];


            for(it_list = (*it_vector).begin(); it_list != (*it_vector).end();++it_list)
            {

                y[LOCATION((*it_list).location.first,(*it_list).location.second,n_pixels)] += phi * (*it_list).probability;

            }

          }
          t++;

        }

        for(int p = 0; p < n_pixels * n_pixels; ++p)
        {

                if(scale[p] > 0){
                rho[p] *= (y[p]/scale[p]) ;
                }
                 //rho[LOCATION(x,p,n_pixels)] = rho[LOCATION(x,p,n_pixels)] *y[LOCATION(x,p,n_pixels)];
        }

    }

    return rho;

  }

  int get_n_pixels(){return n_pixels;}

private:

  int n_tubes;
  int n_pixels;
  int n_iter;
  int emissions;
  int n_lors;
  static std::vector< std::vector<hits_per_pixel> > system_matrix;
  static std::vector<int> list_of_lors;
  static std::vector<int> n;
  float *scale;


  int get_mean_per_lor(file_half &a,file_half &b,std::vector<mean_per_lor> &mean){

    std::vector<mean_per_lor>::iterator it;

    for(it = mean.begin(); it != mean.end(); ++it )
    {

      if(a == (*it).Lor.first && b == (*it).Lor.second)
      {
        printf("Equal: LOR(%d,%d)\n",(*it).Lor.first,(*it).Lor.second);
        return (*it).n;
      }
    }

    return 0.f;
  }

};


