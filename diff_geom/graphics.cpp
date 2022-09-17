#include "graphics.h"


Vec<3> Screen::lookAtPixel(size_t x, size_t y)
{
    double real_x = (double(x) * width) / double(w_resolution) - (0.5 * width);
    double real_y =  (0.5 * height) - (double(y) * height) / double(h_resolution);
    return {real_x, real_y, distance_to_screen}.normalized();
}

Renderer::Renderer(const RiemannianManifold<3> &manifold, genPoint<3> pt) : manifold(manifold), basePoint(pt)
{

}

void Renderer::setBasePoint(genPoint<3> newBasePoint)
{
    basePoint = newBasePoint;
}

void Renderer::setBasis(Vec<3> x_axe_, Vec<3> y_axe_, Vec<3> z_axe_)
{
    auto newBasis = manifold.orthogonalize(basePoint, {x_axe_, y_axe_, z_axe_});
    x_axe = newBasis[0];
    y_axe = newBasis[1];
    z_axe = newBasis[2];
}
