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

    // Projections
    proj_north_plane_to_sphere = [this](Point<N> p){
        Point<N+1> pp;
        for(size_t i = 1; i < N+1; i++) {
            pp[i] = p[i-1];
        }

        pp[0] = 1; // Upper plane!
        Vec<N+1> a_pp = pp - this->south_pole;
        double t = 4 / a_pp.norm2();

        return this->south_pole + a_pp * t;
    };

    proj_south_plane_to_sphere = [this](Point<N> p){
        Point<N+1> pp;
        for(size_t i = 1; i < N+1; i++) {
            pp[i] = p[i-1];
        }

        pp[0] = -1; // Down plane!
        Vec<N+1> a_pp = pp - this->north_pole;
        double t = 4 / a_pp.norm2();

        return this->north_pole + a_pp * t;
    };

    n_proj_to_plane = [](Point<N+1> p1){
        if (p1[0] == -1) return Point<N>::zero();

        double t = 2.0 / (p1[0] + 1.0);

        Point<N> res;
        for(size_t i = 0; i < N; i++) {
            res[i] = t * p1[i+1];
        }
        return res;
    };

    s_proj_to_plane = [](Point<N+1> p1){
        if (p1[0] == 1) return Point<N>::zero();

        double t = 2.0 / (p1[0] - 1.0);

        Point<N> res;
        for(size_t i = 0; i < N; i++) {
            res[i] = t * p1[i+1];
        }
        return res;
    };
    // Structure maps

    this->structureMaps = {
        {{0,1}, [this](Point<N> p){return p[0] == 0.0 ? std::nullopt : std::make_optional( s_proj_to_plane(proj_north_plane_to_sphere(p)) );} },
        {{1,0}, [this](Point<N> p){return p[0] == 0.0 ? std::nullopt : std::make_optional( n_proj_to_plane(proj_south_plane_to_sphere(p)) );} }
    };

    // Metric
    this->metric = {
        InducedMetricTensor<N, N+1>(proj_north_plane_to_sphere),
        InducedMetricTensor<N, N+1>(proj_south_plane_to_sphere)
    };
}

template class Sphere<2>;
template class Sphere<3>;
