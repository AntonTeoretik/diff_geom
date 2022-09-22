#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "algebra.h"
#include "manifold.h"

enum PLAIN_POSITION {NORTH, SOUTH};

template<std::size_t N>
class Inversion : public diffeomorphism<N>
{
protected:
    double radious2;
public:
    Inversion(double radious = 1.0) : radious2(radious*radious) {}
    Point<N> operator ()(const Point<N> &);
};

template<std::size_t N>
class EllipsoidMetric : public InducedMetricTensor<N, N+1>
{
protected:
    double shift;
    std::array<double, N+1> proportions;
public:
    EllipsoidMetric(const std::array<double, N+1>& proportions, const PLAIN_POSITION pos);
    Point<N+1> gen_func(const Point<N> &) const;
};


template<std::size_t N>
class Ellipsoid : public AbstractRiemannianManifold<N>
{
protected:
    double struct_const;
    EllipsoidMetric<N> up_metric, down_metric;
public:
    Ellipsoid(std::array<double, N+1> proportions, double struct_const = 9.0);

    // AbstractManifold interface
public:
    virtual bool changePointIndex(Point<N> &pt, chart_index oldIndex, chart_index newIndex) const;
    virtual bool isPoint(const Point<N> &pt, chart_index i) const;
    const MetricTensor<N>& getMetric(chart_index i) const;
};


#endif // ELLIPSOID_H
