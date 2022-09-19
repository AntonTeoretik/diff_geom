#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "algebra.h"
#include "manifold.h"

template<std::size_t N>
Point<N+1> north_pole = basis<N+1>[0];

template<std::size_t N>
Point<N+1> south_pole = Point<N+1>::zero() - basis<N+1>[0];

///
/// \brief The Sphere class
/// index 0 : North plane
/// index 1 : South plane
///
template<std::size_t N>
class Sphere : public RiemannianManifold<N>
{
public:
    std::function<void(const Point<N+1>&, Point<N>&)> n_proj_to_plane;
    std::function<void(const Point<N+1>&, Point<N>&)> s_proj_to_plane;

    structMap<N> up_to_down, down_to_up;

public:
    Sphere(double controlConst = 3.0);

};

template <size_t N>
void proj_north_plane_to_sphere (const Point<N>&, Point<N+1>&);

template <size_t N>
void proj_south_plane_to_sphere (const Point<N>&, Point<N+1>&);



class Ellipsoid
{
public:
    Ellipsoid();
};

#endif // ELLIPSOID_H
