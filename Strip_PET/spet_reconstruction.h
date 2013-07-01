#ifndef SPET_RECONSTRUCTION_H
#define SPET_RECONSTRUCTION_H

#include <cmath>
#include <vector>
#include <algorithm>
#include "data_structures.h"

//algebra
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

#define INVERSE_PI 0.31830988
#define INVERSE_TWO_PI 0.15915494

using namespace boost::numeric::ublas;



template <typename T = float>
class spet_reconstruction
{
private:
    T R_distance;
    T Scentilator_length;
    std::vector< event<T> > event_list;
    matrix<T> C;
    matrix<T> inverse_C;
    T inverse_C_factor;


public:
    spet_reconstruction(T &R_distance,T &Scentilator_length,T &sigma_z,T &sigma_dl,T &gamma):
        R_distance(R_distance),
        Scentilator_length(Scentilator_length)
    {

        C = matrix<T>(3,3);
        inverse_C = matrix<T>(3,3);

        //matrix C elements

        C(0,0) = sigma_z * sigma_z;
        C(0,1) = static_cast<T>(0);
        C(0,2) = gamma;
        C(1,0) = static_cast<T>(0);
        C(1,1) = sigma_z * sigma_z;
        C(1,2) = - gamma;
        C(2,0) = gamma;
        C(2,1) = -gamma;
        C(2,2) = sigma_dl * sigma_dl;

        //matrix inverse_C elements

        inverse_C(0,0) = (gamma*gamma) - ((sigma_z * sigma_z) * (sigma_dl * sigma_dl));
        inverse_C(0,1) = gamma * gamma;
        inverse_C(0,2) = gamma * (sigma_z * sigma_z);
        inverse_C(1,0) = gamma * gamma;
        inverse_C(1,1) = (gamma*gamma) - ((sigma_z * sigma_z) * (sigma_dl * sigma_dl));
        inverse_C(1,2) = -gamma * (sigma_z * sigma_z);
        inverse_C(2,0) = gamma * (sigma_z * sigma_z);
        inverse_C(2,1) = -gamma * (sigma_z * sigma_z);
        inverse_C(2,2) = -(sigma_z * sigma_z) * (sigma_z * sigma_z);



        T temp_factor = ((gamma*gamma) * (sigma_z*sigma_z)) + ((gamma*gamma) * (sigma_z*sigma_z))+
                        ((sigma_z*sigma_z)*(sigma_z*sigma_z)*(sigma_dl*sigma_dl));
        inverse_C_factor = static_cast<T>(1.0)/temp_factor;

        //axpy_prod(test1, test2, test);

        boost::numeric::ublas::matrix<T> o(1,3);
        boost::numeric::ublas::matrix<T> a(1,3);
        boost::numeric::ublas::matrix<T> b(1,3);
        boost::numeric::ublas::matrix<T> output(1,1);


        matrix<T> test(3,3);
        matrix<T> test2(1,3);

        for(unsigned i = 0; i < test.size1(); ++i){
            for(unsigned j = 0; j < test.size2(); ++j){

                test(i,j) = 1.0;
                //test2(i,j) = 1.0;
            }
        }

        o(0,0) = 1;
        o(0,1) = 1;
        o(0,2) = 1;

        std::cout << o << std::endl;

        a(0,0) = 1;
        a(0,1) = 1;
        a(0,2) = 1;

        b(0,1) = 1;
        b(0,1) = 1;
        b(0,2) = 1;

        axpy_prod(o, test, test2,true);
        std::cout <<"TEST2: " << test2 << std::endl;
        std::cout <<"o: " << trans(o) << std::endl;

        axpy_prod(test2, trans(o), output,true);

        std::cout << output << std::endl;



    };

    T sensitivity(T &y,T &z){

        T a = (0.5 *Scentilator_length - z)/(R_distance - y);
        T b = (0.5 *Scentilator_length + z)/(R_distance + y);
        T c = (0.5 *Scentilator_length + z)/(R_distance - y);

        return INVERSE_PI*(std::atan(std::min(a,b)) - std::atan(std::min(-c,-b)));
    }


    T kernel(int &event,T &y,T &z){

        boost::numeric::ublas::matrix<T> o(1,3);
        boost::numeric::ublas::matrix<T> a(1,3);
        boost::numeric::ublas::matrix<T> b(1,3);
        boost::numeric::ublas::matrix<T> output(1,1);
        boost::numeric::ublas::matrix<T> temp(1,3);

        o(0,0) = 1;
        o(0,1) = 1;
        o(0,2) = 1;

        //std::cout << o << std::endl;

        a(0,0) = 1;
        a(0,1) = 1;
        a(0,2) = 1;

        b(0,0) = 1;
        b(0,1) = 1;
        b(0,2) = 1;


        T vector_element;
        T sqrt_det;
        T exp_full_element;
        T exp_element_3;
        T exp_element_4;

        //człon 1 (a * C^(-1) * a^(T))
        axpy_prod(a, inverse_C, temp,true);
        axpy_prod(temp, trans(a), output,true);

        T temp_a = output(0,0);

        std::cout << output << std::endl;

        //czlon 2 (o * C^(-1) * b^(T))
        axpy_prod(o, inverse_C, temp,true);
        axpy_prod(temp, trans(b), output,true);

        vector_element = temp_a + 2.0*output(0,0);
        std::cout << output << std::endl;
        std::cout << "HERE" << std::endl;
        //człon 3 (b * C^(-1) * a^(T))

        axpy_prod(b, inverse_C, temp,true);
        axpy_prod(temp, trans(a), output,true);

        exp_element_3 = output(0,0);
        std::cout << exp_element_3 << std::endl;


        //czlon 4  człon 3 (b * C^(-1) * b^(T))
        axpy_prod(b, inverse_C, temp,true);
        axpy_prod(temp, trans(b), output,true);

        exp_element_4 = output(0,0);

        std::cout << exp_element_4 << std::endl;

        exp_full_element = exp_element_4 - ((exp_element_3 * exp_element_3)/vector_element);

        //return INVERSE_TWO_PI * (sqrt_det/std::sqrt(vector_element)) * std::exp(-0.5*exp_full_element);
        return 1;
    }


    T tan_angle_from_event(T &z_u,T &z_d) const {return (z_u - z_d)/(2.0*R_distance);}
    T y_from_event(T &dl,T &tan_event) const {return -0.5 * (dl/sqrt(1 + (tan_event*tan_event)));}
    T z_from_event(T &z_u,T &z_d,T &y,T &tan_event) const {return 0.5*(z_u + z_d + (2.0*y *tan_event));}

};

#endif // SPET_RECONSTRUCTION_H
