#include "ellipsoid.h"
#include <iostream>

using namespace std;

#include "algebra.h"
#include "metric.h"
#include "graphics.h"

#include <fstream>

#include <array>

class my_test_funtion_on_manifold : public genFunction<3, Vec<3>>
{
    const Ellipsoid<3>& S;

    double closeness_to_int(double x) const {
        double x1 = std::round(x);
        x1 = abs(x1 - x);
        return x1 <= 0.1 ? 2.0 : 0.0;
    }
    // genFunction interface
public:
    my_test_funtion_on_manifold(const Ellipsoid<3>& S) : S(S) {}

    Vec<3> operator ()(const genPoint<3> & p) const
    {
        Point<3> pp = p.p;

        if(p.i == 1 and !S.changePointIndex(pp, 1, 0)) {
                return {0.0, 0.0, 0.0};
        }
        //pp.scale(0.1);

        double x = closeness_to_int( pp[0] - 0.5);
        double y = closeness_to_int( pp[1] - 0.5);
        double z = closeness_to_int( pp[2] - 0.5);

        return {x*y, y*z, z*x};
    }
};


int main()
{
    Ellipsoid<3> S({1.0, 1.0, 1.0, 1.0});
    my_test_funtion_on_manifold f(S);
    standartWeight w;

    //auto metr = EllipsoidMetric<3>({1,1,1,1}, NORTH);
    //std::cout << "Metr " << metr({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}) << std::endl ;

    Renderer rend(S, f, {0, {0,0,0}}, w);

    rend.number_of_points = 500;
    rend.step = 0.01;
    rend.screen.h_resolution = 50;
    rend.screen.w_resolution = 50;

    rend.setBasis({1,3,2}, {3,2,1}, {1, 0, 0});

    bitmap_image img = rend.render();
    img.save_image("../img.bmp");

    return 0;
}
