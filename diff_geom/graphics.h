#ifndef GRAPHICS_H
#define GRAPHICS_H

#include "manifold.h"
#include "color.h"

struct Screen
{
    double distance_to_screen;
    double width;
    double height;

    size_t w_resolution;
    size_t h_resolution;

    Vec<3> lookAtPixel(size_t x, size_t y) const;
};

class Renderer
{
    const RiemannianManifold<3>& manifold;
    genPoint<3> basePoint;
    Vec<3> x_axe = {1.0, 0.0, 0.0},
           y_axe = {0.0, 1.0, 0.0},
           z_axe = {0.0, 0.0, 1.0};

    Vec<3> trace_one_pixel(const std::function<Vec<3>(genPoint<3>)>& func, size_t x, size_t y) const;

public:
    std::function<double(double)> weight = [](auto){return 1.0;};
    double step = time_step;
    Screen screen = {1.0, 2.0, 2.0, 200, 200};

    size_t number_of_points = 100;
    size_t dist = 1;

    Renderer(const RiemannianManifold<3>& manifold, genPoint<3> basePt);
    void setBasePoint(genPoint<3> newBasePoint);
    void setBasis(Vec<3> x_axe_, Vec<3> y_axe_, Vec<3> z_axe_);

    bitmap_image render(const std::function<Vec<3>(genPoint<3>)>& func);

};

#endif // GRAPHICS_H
