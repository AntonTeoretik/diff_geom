#include "graphics.h"


Vec<3> Screen::lookAtPixel(size_t x, size_t y) const
{
    double real_x = (double(x) * width) / double(w_resolution) - (0.5 * width);
    double real_y =  (0.5 * height) - (double(y) * height) / double(h_resolution);
    return Vec<3>{real_x, real_y, distance_to_screen}.normalized();
}

rgb_t Renderer::trace_one_pixel(const std::function<Vec<3> (genPoint<3>)>& func, size_t x, size_t y) const
{
    Vec<3> v = this->screen.lookAtPixel(x, y);
    auto path = this->manifold.geodesic(basePoint, v, number_of_points, dist, step);
    Vec<3> integr = integrateAlongPath(path, func, weight, step);
    return Color(integr).toRB();
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

bitmap_image Renderer::render(const std::function<Vec<3> (genPoint<3>)> &func)
{
    bitmap_image img(screen.w_resolution, screen.h_resolution);
    for(size_t x = 0; x < screen.w_resolution; x++) {
        for(size_t y = 0; y < screen.h_resolution; y++) {
            std::cout << x << " " << y << std::endl;
            rgb_t res = trace_one_pixel(func, x, y);
            std::cout << int(res.red) << " " << int(res.green) << " " << int(res.blue) <<std::endl;
            img.set_pixel(x, y, trace_one_pixel(func, x, y));
        }
    }
    return img;
}
