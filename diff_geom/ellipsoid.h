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
    CUDA_F Inversion(double radious = 1.0) : radious2(radious*radious) {}
    CUDA_F Point<N> operator ()(const Point<N> &);
};

template<std::size_t N>
class EllipsoidMetric : public InducedMetricTensor<N, N+1>
{
protected:
    double shift;
    double proportions[N+1];
public:
    CUDA_F EllipsoidMetric(const Vec<N+1>& proportions, const PLAIN_POSITION pos);
    CUDA_F Point<N+1> gen_func(const Point<N> &) const override;
};


template<std::size_t N>
class Ellipsoid : public AbstractRiemannianManifold<N>
{
protected:
    double struct_const;
    EllipsoidMetric<N> up_metric, down_metric;
public:
    CUDA_F Ellipsoid(const Vec<N+1>& proportions, double struct_const = 9.0);

    // AbstractManifold interface
public:
    virtual CUDA_F bool changePointIndex(Point<N> &pt, chart_index oldIndex, chart_index newIndex) const;
    virtual CUDA_F bool isPoint(const Point<N> &pt, chart_index i) const;
    CUDA_F const MetricTensor<N>& getMetric(chart_index i) const;
};


#endif // ELLIPSOID_H
