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

    Sphere<1> S;


    //auto met = InducedMetricTensor<2,3>(S.proj_north_plane_to_sphere, 0.001);
    auto met = S.metric[0]; // Doesn't work

    //auto m = InducedMetricTensor<2,3>([](Point<2> p){return Point<3>{1, p[0], p[1]};});

    //std::shared_ptr<MetricTensor<2>> met = std::make_shared<InducedMetricTensor<2,3>>(m);

    //std::cout << (*met)({0,0}, {1.0,0.0}, {1.0, 0.0}) << std::endl;

    auto vec = S.geodesic({0, {0}}, {1}, 10000, 0.01);
    for (size_t i = 0; i < vec.size(); i++) {
        if (i % 10 == 0) {
            auto pt = vec[i];
            std::cout << pt.i << ": " << pt.p.to_str() << std::endl;
        }
    }


    return 0;
}
