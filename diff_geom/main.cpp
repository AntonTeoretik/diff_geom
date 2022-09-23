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
        /*if (pp.norm2() >= 100) {
            return {100.0, 100.0, 100.0};
        }*/
        pp.scale(4);

        int res = int(0.1 - abs(pp[0]) > 0) + int(0.1 - abs(pp[1]) > 0) + int(0.1 - abs(pp[2]) > 0);
        if (res < 1) {
            return {0.0, 0.0, 0.0};
        }
        double x = closeness_to_int(pp[0]);
        double y = closeness_to_int(pp[1]);
        double z = closeness_to_int(pp[2]);
        //double z = std::max(0.0, 0.1 - abs(pp[2]));

        return Vec<3>{x*y, y*z, z*x}*3.6;
    }
};


int main()
{
    Ellipsoid<3> S({1.0, 1.0, 1.0, 1.0});
    my_test_funtion_on_manifold f(S);
    standartWeight w(0.2);

    //auto metr = EllipsoidMetric<3>({1,1,1,1}, NORTH);
    //std::cout << "Metr " << metr({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}) << std::endl ;

    Renderer rend(S, f, {0, {5,5,1}}, w);


    rend.number_of_points = 330;
    rend.step = 0.01;
    rend.screen.h_resolution = 9  * 12;
    rend.screen.w_resolution = 16 * 12;
    //1366x768

    rend.screen.height = 7.68;
    rend.screen.width = 13.66;
    rend.screen.distance_to_screen = 3.0;


    rend.setBasis({-5, -5, -1}, {0, 0, -1}, {0, 1, 0});



    bitmap_image img = rend.render();
    img.save_image("../img.bmp");

    return 0;
}
