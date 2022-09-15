#include <iostream>

using namespace std;

#include "algebra.h"
#include "metric.h"
#include "manifold.h"
#include "ellipsoid.h"

#include <fstream>

#include <array>

int main()
{

    Sphere<2> S;


    auto met = InducedMetricTensor<2,3>(S.proj_north_plane_to_sphere, 0.001);
    //auto met = S.metric[0]; // Doesn't work

    std::cout << met({0,0}, {1.0,0.0}, {1.0, 0.0}) << std::endl;



    return 0;
}
