#include "graphics.h"


Vec<3> Screen::lookAtPixel(size_t x, size_t y) const
{
    double real_x = (double(x) * width) / double(w_resolution) - (0.5 * width);
    double real_y =  (0.5 * height) - (double(y) * height) / double(h_resolution);
    return Vec<3>{real_x, real_y, distance_to_screen}.normalized();
}

Vec<3> Renderer::trace_one_pixel(const std::function<Vec<3> (genPoint<3>)>& func, size_t x, size_t y) const
{
    Vec<3> v = this->screen.lookAtPixel(x, y);
    Vec<3> v1 = x_axe * v[0] + y_axe * v[1] + z_axe * v[2];
    return this->manifold.integrateAlongPath(basePoint, v1, number_of_points,
                                                 func, weight, step);
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
    auto newBasis = manifold.orthogonalize(basePoint, {x_axe_, y_axe_, z_axe_}, true);
    x_axe = newBasis[0];
    y_axe = newBasis[1];
    z_axe = newBasis[2];

    std::cout << "Basis set:" << std::endl;
    std::cout << x_axe.to_str() << std::endl;
    std::cout << y_axe.to_str() << std::endl;
    std::cout << z_axe.to_str() << std::endl;
}

bitmap_image Renderer::render(const std::function<Vec<3> (genPoint<3>)> &func)
{
    bitmap_image img(screen.w_resolution, screen.h_resolution);

    std::vector<std::vector<Vec<3>>> table_of_integrals(screen.w_resolution);

    for(size_t x = 0; x < screen.w_resolution; x++) {
        for(size_t y = 0; y < screen.h_resolution; y++) {
            Vec<3> res = trace_one_pixel(func, x, y);
            table_of_integrals[x].push_back(res);
        }
        std::cout << x << std::endl;
    }

    for(size_t x = 0; x < screen.w_resolution; x++) {
        for(size_t y = 0; y < screen.h_resolution; y++) {
            img.set_pixel(x, y, Color(table_of_integrals[x][y]).toRGB());
        }
    }


    return img;
}
