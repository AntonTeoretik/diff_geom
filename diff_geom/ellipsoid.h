#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "algebra.h"
#include "manifold.h"


///
/// \brief The Sphere class
/// index 0 : North plane
/// index 1 : South plane
///
template<std::size_t N>
class Sphere : public RiemannianManifold<N>
{
public:
    Point<N+1> north_pole;
    Point<N+1> south_pole;

    std::function<void(const Point<N>&, Point<N+1>&)> proj_north_plane_to_sphere;
    std::function<void(const Point<N>&, Point<N+1>&)> proj_south_plane_to_sphere;

    std::function<void(const Point<N+1>&, Point<N>&)> n_proj_to_plane;
    std::function<void(const Point<N+1>&, Point<N>&)> s_proj_to_plane;

    structMap<N> up_to_down, down_to_up;

public:
    Sphere(double controlConst = 3.0);

    Point<N+1> get_north_pole() const;
    Point<N+1> get_south_pole() const;
};




class Ellipsoid
{
public:
    Ellipsoid();
};

#endif // ELLIPSOID_H
