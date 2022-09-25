#ifndef GRAPHICS_H
#define GRAPHICS_H

#include "manifold.h"
#include "color.h"

#include <omp.h>
#include "util.h"

struct Screen
{
    double distance_to_screen;
    double width;
    double height;

    size_t w_resolution;
    size_t h_resolution;

    CUDA_F Vec<3> lookAtPixel(size_t x, size_t y) const;
};

class Renderer
{
    genPoint<3> basePoint;
    Vec<3> x_axe;
    Vec<3> y_axe;
    Vec<3> z_axe;


public:
    CUDA_F Vec<3> trace_one_pixel(size_t x, size_t y) const;

    const weightFunction& weight;
    const genFunction<3, Vec<3>>& func;
    const AbstractRiemannianManifold<3>& manifold;

    double step = time_step;
    Screen screen = {1.0, 2.0, 2.0, 200, 200};
    size_t number_of_points = 100;

    CUDA_F Renderer(const AbstractRiemannianManifold<3>& manifold,
             const genFunction<3, Vec<3>>& func,
             const genPoint<3>& basePt,
             const weightFunction& weight = standartWeight());

    CUDA_F void setBasePoint(genPoint<3> newBasePoint);
    CUDA_F void setBasis(Vec<3> x_axe_, Vec<3> y_axe_, Vec<3> z_axe_);

    bitmap_image render();
};

#endif // GRAPHICS_H
