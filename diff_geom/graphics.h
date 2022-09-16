#ifndef GRAPHICS_H
#define GRAPHICS_H

#include "manifold.h"


struct Screen
{
    double distance_to_screen;
    double width;
    double height;

    size_t w_resolution;
    size_t h_resolution;

    Vec<3> lookAtPixel(size_t x, size_t y);
};

class Renderer
{
    const RiemannianManifold<3>& manifold;
    genPoint<3> basePoint;
    Vec<3> x_axe = {1.0, 0.0, 0.0},
           y_axe = {0.0, 1.0, 0.0},
           z_axe = {0.0, 0.0, 1.0};
public:

    Screen screen = {1.0, 2.0, 2.0, 400, 400};

    Renderer(const RiemannianManifold<3>& manifold, genPoint<3> basePt);
    void setBasePoint(genPoint<3> newBasePoint);
    void setBasis(Vec<3> x_axe_, Vec<3> y_axe_, Vec<3> z_axe_);
};

#endif // GRAPHICS_H
