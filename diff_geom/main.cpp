#include <iostream>

using namespace std;

#include <algebra.h>
#include <geometry.h>

#include <fstream>

#include <array>

int main()
{
    MetricTensor<2> g([](Point<2> P, Vec<2> v1, Vec<2> v2){ // polar coordinates
        Matrix2D<2> M = { {1, 0 },
                          {0, P[0]*P[0]} };
        return v1 * (M * v2);
    });


    MetricTensor<2> h;

    std::function<Point<3>(Point<2>)> surface = [](Point<2> p){ return Point<3>{p[0], p[1], p[0]*p[0] + p[1]*p[1]}; };
    InducedMetricTensor surf(surface);

    auto geod = geodesic(surf, {2, 0}, {0, 1}, 100000);

    auto fout = std::ofstream("../test.csv");

    int i = 0;
    fout << "x,y,z\n" ;
    for (auto p : geod) {
        i++;
        if(i % 100 == 0) {
            fout << p[0] << "," << p[1] << "," << p[0]*p[0] + p[1]*p[1] << "\n";
            i = 1;
        }
    }
    fout.close();

    return 0;
}
