#include "ellipsoid.h"


template<std::size_t N>
void proj_north_plane_to_sphere(const Point<N> & p, Point<N+1> & pp)
{
    for(size_t i = 1; i < N+1; i++) {
        pp[i] = p[i-1];
    }

    pp[0] = 1; // Upper plane!
    Vec<N+1> a_pp = pp;
    a_pp[0] += 1;

    double t = 4.0 / a_pp.norm2();
    a_pp.scale(t);
    pp = a_pp;
    pp[0] -= 1;
}

template<std::size_t N>
void proj_south_plane_to_sphere(const Point<N> & p, Point<N+1> & pp)
{
    for(size_t i = 1; i < N+1; i++) {
        pp[i] = p[i-1];
    }
    pp[0] = -1; // Down plane!

    Vec<N+1> a_pp = pp;

    a_pp[0] -= 1;
    double t = 4.0 / a_pp.norm2();
    a_pp.scale(t);
    pp = a_pp;
    pp[0] += 1;
}

template<std::size_t N>
Sphere<N>::Sphere(double controlConst) :
    RiemannianManifold<N>({}, {}, {})
{
    // Atlas
    Chart<N> open_circle = [controlConst](Point<N> p){return (p.norm() < controlConst);};
    this->atlas = {open_circle, open_circle};
    this->atlas_size = 2;


    // Structure maps
    this->structureMaps = {
        {{0,1},
         [](Point<N> p){
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
         [](Point<N> p){

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
        std::make_shared<InducedMetricTensor<N, N+1>>(proj_north_plane_to_sphere<N>),
        std::make_shared<InducedMetricTensor<N, N+1>>(proj_south_plane_to_sphere<N>)
    };
}

template<size_t N>
void s_proj_to_plane(const Point<N+1> & p1, Point<N> & res)
{
    if (p1[0] == 1) {
        res = Point<N>::zero();
        return;
    }
    double t = -2.0 / (p1[0] - 1.0);
    for(size_t i = 0; i < N; i++) {
        res[i] = t * p1[i+1];
    }
}

template<size_t N>
void n_proj_to_plane(const Point<N+1> & p1, Point<N> & res)
{
    if (p1[0] == -1) {
        res = Point<N>::zero();
        return;
    }

    double t = 2.0 / (p1[0] + 1.0);

    for(size_t i = 0; i < N; i++) {
        res[i] = t * p1[i+1];
    }
}



template class Sphere<1>;
template class Sphere<2>;
template class Sphere<3>;


