#include <iostream>

using namespace std;

#include "algebra.h"
#include "metric.h"
#include "manifold.h"
#include "ellipsoid.h"
#include "graphics.h"

#include <fstream>

#include <array>


double closeness_to_int(double x) {
    double x1 = std::round(x);
    return 2*(0.5 - abs(x1 - x));
}

Vec<3> my_func(const genPoint<3>& p, const Sphere<3>& S) {
    Point<3> pp = p.p;

    if(p.i == 1) {
        auto res = S.changePointIndex(p, 0);
        if (!res.has_value())
            return {0.0, 0.0, 0.0};
        pp = res.value();
    }
    pp.scale(0.1);

    double x = closeness_to_int( pp[0] - 0.5);
    double y = closeness_to_int( pp[1] - 0.5);
    double z = closeness_to_int( pp[2] - 0.5);

    return {x*y, y*z, z*x};

}

int main()
{
    double ts = 0.01;
    Sphere<3> S(3.0);

    Renderer rend(S, {1, {0,0,0}});

    rend.number_of_points = 100;
    rend.step = 0.05;
    rend.screen.h_resolution = 200;
    rend.screen.w_resolution = 200;

    rend.weight = [](double x){return 1 / (x + 1);};

    bitmap_image img = rend.render([&S](auto p){return my_func(p, S);});
    img.save_image("../img.bmp");

    return 0;
}
