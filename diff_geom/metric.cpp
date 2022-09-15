#include "metric.h"



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
double MetricTensor<N>::operator()(Point<N> p, Vec<N> vec1, Vec<N> vec2) const
{
    return g(p, vec1, vec2);
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
Matrix2D<N> MetricTensor<N>::kristMatrix(std::size_t l, Point<N> p) const
{
    std::cout << "Krist? " << p.to_str() << std::endl;
    return Matrix2D<N>([this, l, p](auto i, auto j){return this->krist(l, i, j, p);});
}


template class MetricTensor<1>;
template class MetricTensor<2>;
template class MetricTensor<3>;


template<std::size_t N>
std::vector<Point<N> > geodesic(const MetricTensor<N> &g, Point<N> start, Vec<N> vel, std::size_t steps_num, double dt)
{
    Point<N> prev = start;
    Point<N> now = start + (vel * dt);
    Point<N> next;

    std::vector<Point<N>> res = {prev, now};

    for (std::size_t i = 0; i < steps_num; i++) {
        res.push_back(now);

        Vec<N> curv;
        for (size_t l = 0; l < N; l++)
            curv[l] = (now - prev) * (g.kristMatrix(l, now) * (now - prev));

        next = now * 2 - prev - curv;

        prev = now;
        now = next;
    }

    return res;
}



template std::vector<Point<2> > geodesic(const MetricTensor<2> &g, Point<2> start, Vec<2> vel, std::size_t steps_num, double dt);
template std::vector<Point<3> > geodesic(const MetricTensor<3> &g, Point<3> start, Vec<3> vel, std::size_t steps_num, double dt);


template<std::size_t N, std::size_t M>
InducedMetricTensor<N, M>::InducedMetricTensor(std::function<Point<M> (Point<N>)> f, double pres) : f(f), pres(pres)
{
    // This is ALMOST linear version
    this->g = [this](Point<N> p, Vec<N> v1, Vec<N> v2)
    {
        double alpha1 = v1.norm();
        double alpha2 = v2.norm();

        if(alpha1 * alpha2 == 0) return 0.0; // Almost impossible, but anyway

        std::cout << v1.to_str() << "; " << v2.to_str() << std::endl;

        Vec<M> dv1 = this->f(p + v1.nomalized() * this->pres) - this->f(p - v1.nomalized() * this->pres);
        Vec<M> dv2 = this->f(p + v2.nomalized() * this->pres) - this->f(p - v2.nomalized() * this->pres);

        std::cout << dv1.to_str() << "; " << dv2.to_str() << "; " << this->pres << std::endl;
        std::cout << alpha1 << "; " << alpha2 << ";" << dv1 * dv2 << "; " << this->pres << std::endl;
        return 0.25 * alpha1 * alpha2 * (dv1 * dv2) / (this->pres * this->pres);
    };
}

template<std::size_t N, std::size_t M>
Point<M> InducedMetricTensor<N, M>::apply_generator(Point<N> p) const
{
    return f(p);
}

template class InducedMetricTensor<1,2>;
template class InducedMetricTensor<2,3>;
template class InducedMetricTensor<3,4>;
template class InducedMetricTensor<2,2>;
