#include "ellipsoid.h"

template<std::size_t N>
CUDA_F Point<N> Inversion<N>::operator ()(const Point<N> & P)
{
    double norm2 = P.norm2();
    if(norm2 == 0) {
        return P;
    }
    Point<N> P1 = P;
    P1.scale(radious2 / norm2);
    return P1;
}

template<std::size_t N>
CUDA_F EllipsoidMetric<N>::EllipsoidMetric(const Vec<N+1>& _proportions, const PLAIN_POSITION pos) 
{
    for (int i = 0; i < N + 1; ++i) {
        this->proportions[i] = _proportions[i];
    }
    if(pos == NORTH) {
        shift = -1.0;
    } else if (pos == SOUTH){
        shift = 1.0;
    }
}

template<std::size_t N>
CUDA_F Point<N+1> EllipsoidMetric<N>::gen_func(const Point<N> & p) const
{
    Point<N+1> pp;
    for(size_t i = 1; i < N+1; i++) {
        pp[i] = p[i-1];
    }

    pp[0] = shift;
    Vec<N+1> a_pp = pp;
    a_pp[0] += shift;

    double t = 4.0 / a_pp.norm2();
    a_pp.scale(t);
    pp = a_pp;
    pp[0] -= shift;

    for(size_t i = 0; i < N+1; i++) {
        pp[i] *= proportions[i];
    }

    return pp;
}

template<std::size_t N>
CUDA_F Ellipsoid<N>::Ellipsoid(const Vec<N+1>& proportions, double struct_const) :
    struct_const(struct_const),
    up_metric(EllipsoidMetric<N>(proportions, NORTH)),
    down_metric(EllipsoidMetric<N>(proportions, SOUTH))
{
    this->atlas_size = 2;
}

template<std::size_t N>
CUDA_F bool Ellipsoid<N>::changePointIndex(Point<N> &pt, chart_index oldIndex, chart_index newIndex) const
{
    //std::cout << "Ellipsoid<N>::changePointIndex: " << pt.to_str() << "; " << oldIndex << " " << newIndex << std::endl;
    if (oldIndex == newIndex) {
        return false;
    }

    double norm2 = pt.norm2();
    if(norm2 == 0.0) {
        return false;
    }
    pt.scale(4.0 / norm2);
    return true;
}

template<std::size_t N>
CUDA_F bool Ellipsoid<N>::isPoint(const Point<N> &pt, chart_index) const
{
    return pt.norm2() <= struct_const;
}

template<std::size_t N>
CUDA_F const MetricTensor<N> &Ellipsoid<N>::getMetric(chart_index i) const
{
    if (i == 0) return this->up_metric;
    if (i == 1) return this->down_metric;
    return this->down_metric;

    // exit(1);
    // throw std::invalid_argument("Ellipsoid<N>::getMetric : chart index must be either 0 or 1");
}

template class Ellipsoid<3>;
template class EllipsoidMetric<3>;
