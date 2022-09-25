#include "graphics.h"


CUDA_F Vec<3> Screen::lookAtPixel(size_t x, size_t y) const
{
    double real_x = (double(x) * width) / double(w_resolution) - (0.5 * width);
    double real_y =  (0.5 * height) - (double(y) * height) / double(h_resolution);

    // return Vec<3>{real_x, real_y, distance_to_screen}.normalized();

    return Vec<3>{distance_to_screen, real_x, real_y};
}

CUDA_F Vec<3> Renderer::trace_one_pixel(size_t x, size_t y) const
{
    Vec<3> v = this->screen.lookAtPixel(x, y);
    Vec<3> v1 = x_axe * v[0] + y_axe * v[1] + z_axe * v[2];

    return manifold.integrateAlongPath(basePoint,
                                       v1,
                                       number_of_points,
                                       func,
                                       weight,
                                       step);
}

CUDA_F Renderer::Renderer(const AbstractRiemannianManifold<3> &manifold,
                   const genFunction<3, Vec<3> > &func,
                   const genPoint<3> &basePt,
                   const weightFunction &weight) :
    basePoint(basePt),
    weight(weight),
    func(func),
    manifold(manifold)
{
    x_axe = Vec<3>::zero();
    x_axe[0] = 1.0;
    y_axe = Vec<3>::zero();
    y_axe[1] = 1.0;
    z_axe = Vec<3>::zero();
    z_axe[2] = 1.0;
}


CUDA_F void Renderer::setBasePoint(genPoint<3> newBasePoint)
{
    basePoint = newBasePoint;
}

CUDA_F void Renderer::setBasis(Vec<3> x_axe_, Vec<3> y_axe_, Vec<3> z_axe_)
{
    Vec<3> vecs[3];
    vecs[0] = x_axe_;
    vecs[1] = y_axe_;
    vecs[2] = z_axe_;
    
    Vec<3> vecs_out[3];
    Span<Vec<3>> in_span{vecs, 3};
    Span<Vec<3>> out_span{vecs_out, 3};
    int n_out = 0;

    manifold.orthogonalize(basePoint, in_span, true, out_span, n_out);

    x_axe = vecs_out[0];
    y_axe = vecs_out[1];
    z_axe = vecs_out[2];

    // std::cout << "Basis set:" << std::endl;
    // std::cout << x_axe.to_str() << std::endl;
    // std::cout << y_axe.to_str() << std::endl;
    // std::cout << z_axe.to_str() << std::endl << std::endl;
}


bitmap_image Renderer::render()
{
    bitmap_image img(screen.w_resolution, screen.h_resolution);

    std::vector<std::vector<Vec<3>>> table_of_integrals(screen.w_resolution);
    for(size_t x = 0; x < screen.w_resolution; x++) {
        table_of_integrals[x].resize(screen.h_resolution);
    }

    for(size_t x = 0; x < screen.w_resolution; x++) {
        for(size_t y = 0; y < screen.h_resolution; y++) {
            Vec<3> res = trace_one_pixel(x, y);
            table_of_integrals[x][y] = res;
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
