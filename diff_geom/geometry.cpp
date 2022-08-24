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
Matrix2D<N> MetricTensor<N>::getMatrix(Point<N> P) const
{
    return Matrix2D<N>([this, P](int i, int j){return this->g(P, basis<N>[i], basis<N>[j]); });
}

template<std::size_t N>
double MetricTensor<N>::getCoord(Point<N> P, std::size_t i, std::size_t j) const
{
    return g(P, basis<N>[i], basis<N>[j]);
}

template<std::size_t N>
double MetricTensor<N>::krist(std::size_t l, std::size_t j, std::size_t k, Point<N> p) const
{
    double res = 0;

    Matrix2D<N> Mg = getMatrix(p);
    Matrix2D<N> MgInv = Mg.inverse();

    for(std::size_t r = 0; r < N; r++) {
        double dk_grj = partialDer<N>(p, basis<N>[k], [this, r, j](Point<N> pp){return this->getCoord(pp, r, j);} );
        double dj_grk = partialDer<N>(p, basis<N>[j], [this, r, k](Point<N> pp){return this->getCoord(pp, r, k);} );
        double dr_gjk = partialDer<N>(p, basis<N>[r], [this, k, j](Point<N> pp){return this->getCoord(pp, j, k);} );

        res += MgInv[l][r] * ( dk_grj + dj_grk - dr_gjk );
    }
    return res * 0.5;
}


template<std::size_t N>
double MetricTensor<N>::operator()(Point<N> p, Vec<N> vec1, Vec<N> vec2) const
{
    return g(p, vec1, vec2);
}

template class MetricTensor<3>;
template class MetricTensor<2>;

