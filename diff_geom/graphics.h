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
    genPoint<3> basePoint;
    Vec<3> x_axe = {1.0, 0.0, 0.0},
           y_axe = {0.0, 1.0, 0.0},
           z_axe = {0.0, 0.0, 1.0};

    Vec<3> trace_one_pixel(size_t x, size_t y) const;

public:
    const weightFunction& weight;
    const genFunction<3, Vec<3>>& func;
    const AbstractRiemannianManifold<3>& manifold;

    double step = time_step;
    Screen screen = {1.0, 2.0, 2.0, 200, 200};
    size_t number_of_points = 100;

    Renderer(const AbstractRiemannianManifold<3>& manifold,
             const genFunction<3, Vec<3>>& func,
             const genPoint<3>& basePt,
             const weightFunction& weight = standartWeight());

    void setBasePoint(genPoint<3> newBasePoint);
    void setBasis(Vec<3> x_axe_, Vec<3> y_axe_, Vec<3> z_axe_);

    bitmap_image render();

};

#endif // GRAPHICS_H
