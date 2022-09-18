#include <iostream>

using namespace std;

#include "algebra.h"
#include "metric.h"
#include "manifold.h"
#include "ellipsoid.h"
#include "graphics.h"

#include <fstream>

#include <array>

int main()
{
    double ts = 0.01;
    Sphere<3> S(3.0);

    Renderer rend(S, {0, {0,0,0}});

    rend.number_of_points = 100;
    rend.step = 0.01;
    rend.screen.h_resolution = 100;
    rend.screen.w_resolution = 100;

    bitmap_image img = rend.render([](auto){return Vec<3>{1.0,1.0,1.0};});
    img.save_image("../img.bmp");

    return 0;
}
