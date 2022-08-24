#include "geometry.h"



template<std::size_t N>
MetricTensor<N>::MetricTensor() : g([](Point<N>, Vec<N> v1, Vec<N> v2) {
    double res = 0;
    for(size_t i = 0; i < N; i++) {
        res += v1[i] * v2[i];
    }
    return res;
})
{

}

template<std::size_t N>
MetricTensor<N>::MetricTensor(std::function<double (Point<N>, Vec<N>, Vec<N>)> g_) : g(g_)
{

}

template<std::size_t N>
double MetricTensor<N>::operator()(Point<N> p, std::array<Vec<N>, 2> vecs) const
{
    return g(p, vecs[0], vecs[1]);
}

template<std::size_t N>
double MetricTensor<N>::krist(std::size_t l, std::size_t j, std::size_t k, Point<N> p) const
{
    double res = 0;


    return res;
}


template<std::size_t N>
double MetricTensor<N>::operator()(Point<N> p, Vec<N> vec1, Vec<N> vec2) const
{
    return g(p, vec1, vec2);
}

template class MetricTensor<3>;
template class MetricTensor<2>;

