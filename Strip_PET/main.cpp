#include <iostream>
#include <vector>

#include "strip_pet.h"
#include "spet_reconstruction.h"
#include "data_structures.h"

using namespace std;

int main()
{
    float R_distance = 10.0f;
    float Scentilator_length = 10.0f;
    float x = 0.0f;
    float y = 0.0f;
    float a = 4.0f;
    float b = 3.0f;
    float phi = 0.0f;

    Strip_PET<> test(R_distance,Scentilator_length,x,y,a,b,phi);

    test.emit_event();

    float sigma = 1;
    float dl = 1;
    float gamma = 0;

    spet_reconstruction<> reconstruction(R_distance,Scentilator_length,sigma,dl,gamma);

    int bleh = 1;
    reconstruction.kernel(bleh,y,a);

    return 0;
}

