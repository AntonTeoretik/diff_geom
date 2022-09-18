#include "ellipsoid.h"

template<std::size_t N>
Point<N+1> Sphere<N>::get_north_pole() const
{
    return north_pole;
}

template<std::size_t N>
Point<N+1> Sphere<N>::get_south_pole() const
{
    return south_pole;
}

template<std::size_t N>
Sphere<N>::Sphere(double controlConst) :
    RiemannianManifold<N>({}, {}, {})
{
    north_pole = Point<N+1>::zero(); north_pole[0] = 1;
    south_pole = Point<N+1>::zero(); south_pole[0] = -1;

    // Atlas
    Chart<N> open_circle = [controlConst](Point<N> p){return (p.norm() < controlConst);};
    this->atlas = {open_circle, open_circle};
    this->atlas_size = 2;

    // Projections
    proj_north_plane_to_sphere = [this](const Point<N>& p, Point<N+1>& pp){
        for(size_t i = 1; i < N+1; i++) {
            pp[i] = p[i-1];
        }

        pp[0] = 1; // Upper plane!
        Vec<N+1> a_pp = pp - this->south_pole;
        double t = 4 / a_pp.norm2();

        pp = (this->south_pole + (a_pp * t));
    };

    proj_south_plane_to_sphere = [this](const Point<N>& p, Point<N+1>& pp){
        for(size_t i = 1; i < N+1; i++) {
            pp[i] = p[i-1];
        }

        pp[0] = -1; // Down plane!
        Vec<N+1> a_pp = pp - this->north_pole;
        double t = 4 / a_pp.norm2();

        pp = this->north_pole + (a_pp * t);
    };

    n_proj_to_plane = [](const Point<N+1>& p1, Point<N>& res){
        if (p1[0] == -1) {
            res = Point<N>::zero();
            return;
        }

        double t = 2.0 / (p1[0] + 1.0);

        for(size_t i = 0; i < N; i++) {
            res[i] = t * p1[i+1];
        }
    };

    s_proj_to_plane = [](const Point<N+1>& p1, Point<N>& res){
        if (p1[0] == 1) {
            res = Point<N>::zero();
            return;
        }
        double t = -2.0 / (p1[0] - 1.0);
        for(size_t i = 0; i < N; i++) {
            res[i] = t * p1[i+1];
        }
    };

    // Structure maps
    this->structureMaps = {
        {{0,1},
         [this](Point<N> p){
             if (p == Point<N>::zero()) {
                 return std::optional<Point<N>>{};
             }
             Point<N+1> pp;
             proj_north_plane_to_sphere(p, pp);
             s_proj_to_plane(pp, p);
             return std::make_optional<Point<N>>(p);
         }
        },

        {{1,0},
         [this](Point<N> p){

             if (p == Point<N>::zero()) {
                 return std::optional<Point<N>>{};
             }
             Point<N+1> pp;
             proj_south_plane_to_sphere(p, pp);
             n_proj_to_plane(pp, p);
             return std::make_optional<Point<N>>(p);
            }
         }
    };

    // Metric
    this->metric = {
        std::make_shared<InducedMetricTensor<N, N+1>>(this->proj_north_plane_to_sphere),
        std::make_shared<InducedMetricTensor<N, N+1>>(this->proj_south_plane_to_sphere)
    };
}



template class Sphere<1>;
template class Sphere<2>;
template class Sphere<3>;
