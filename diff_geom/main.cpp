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
    auto vec = S.geodesic({0, {0}}, {1}, 1000, 0.01);
    /* for (size_t i = 0; i < vec.size(); i++) {
        if (i % 1 == 0) {
            auto pt = vec[i];
            std::cout << pt.i << ": " << pt.p.to_str() << std::endl;
        }
    } */


    return 0;
}
