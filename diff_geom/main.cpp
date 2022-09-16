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

    Sphere<3> S(3.0);
    auto vec = S.geodesic({0, {0,0,0}}, {0,0.1,1}, 1000, 1, 0.01);


    std::ofstream fout("../test.csv");
    fout << "x,y,z" << std::endl;

    //std::cout << S.proj_north_plane_to_sphere({1.4142,1.4142}).to_str() << std::endl;

    for (auto pt : vec) {
        auto new_pt = S.immerse(pt).value();
        //std::cout << pt.i << ": " << pt.p.to_str() << ": " << new_pt.to_str() << std::endl;
        fout << new_pt[0] << "," << new_pt[1] << "," << new_pt[2] << std::endl;
        //fout << pt.p[0] << "," << pt.p[1] << "," << 0 << std::endl;
    }

    fout.close();



    return 0;
}
