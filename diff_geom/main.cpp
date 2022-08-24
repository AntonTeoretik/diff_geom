#include <iostream>

using namespace std;

#include <algebra.h>
#include <geometry.h>

#include <array>

int main()
{
    MetricTensor<2> g([](Point<2> P, Vec<2> v1, Vec<2> v2){ // polar coordinates
        Matrix2D<2> M = { {1, 0 },
                          {0, P[0]*P[0]} };
        return v1 * (M * v2);
    });


    MetricTensor<2> h;

    auto geod = geodesic(g, {1,3.1415926 / 2}, {0, -1}, 100000);
    for (auto p : geod) {
        std::cout << p.to_str() << std::endl;
    }

    return 0;
}
