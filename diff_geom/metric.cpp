#include "metric.h"



template<std::size_t N>
MetricTensor<N>::MetricTensor() : g([](Point<N>&, const Vec<N>& v1, const Vec<N>& v2) {
    double res = 0;
    for(size_t i = 0; i < N; i++) {
        res += v1[i] * v2[i];
    }
    return res;
})
{

}

template<std::size_t N>
MetricTensor<N>::MetricTensor(std::function<double (Point<N>&, const Vec<N>&, const  Vec<N>&)> g_) : g(g_)
{

}

template<std::size_t N>
double MetricTensor<N>::operator()(Point<N>& p, std::array<Vec<N>, 2>& vecs) const
{
    return g(p, vecs[0], vecs[1]);
}


template<std::size_t N>
double MetricTensor<N>::operator()(Point<N>& p, Vec<N>& vec1, Vec<N>& vec2) const
{
    return g(p, vec1, vec2);
}


template<std::size_t N>
Matrix2D<N> MetricTensor<N>::getMatrix(Point<N>& P)
{
    return Matrix2D<N>([this, &P](int i, int j){return this->getCoord(P, i, j); });
}

template<std::size_t N>
double MetricTensor<N>::getCoord(Point<N>& P, std::size_t i, std::size_t j)
{
    return g(P, basis<N>[i], basis<N>[j]);
}

template<std::size_t N>
double MetricTensor<N>::krist(std::size_t l, std::size_t j, std::size_t k, Point<N>& p)
{
    double res = 0;

    Matrix2D<N> MgInv = getMatrix(p).inverse();

    for(std::size_t r = 0; r < N; r++) {
        double dk_grj = partialDer<N>(p, k, [this, r, j](Point<N>& pp){return this->getCoord(pp, r, j);} );
        double dj_grk = partialDer<N>(p, j, [this, r, k](Point<N>& pp){return this->getCoord(pp, r, k);} );
        double dr_gjk = partialDer<N>(p, r, [this, k, j](Point<N>& pp){return this->getCoord(pp, j, k);} );

        res += MgInv[l][r] * ( dk_grj + dj_grk - dr_gjk );
    }
    return res * 0.5;
}

template<std::size_t N>
Matrix2D<N> MetricTensor<N>::kristMatrix(std::size_t l, Point<N>& p)
{
    Matrix2D<N> res;
    for (size_t i = 1; i < N; i++) {
        for (size_t j = 0; j < i; j++) {
            double kr = krist(l, i, j, p);
            res[i][j] = kr;
            res[j][i] = kr;
        }
    }
    for (size_t i = 0; i < N; i++) {
        res[i][i] = this->krist(l, i, i, p);
    }

    //return Matrix2D<N>([this, l, &p](auto i, auto j){return this->krist(l, i, j, p);});
    return res;
}

template<std::size_t N>
std::vector<Vec<N> > MetricTensor<N>::orthogonalize(Point<N>& pt, const std::vector<Vec<N> > &vecs, bool norm) const
{
    std::vector<Vec<N> > res = {};
    for(size_t i = 0; i < vecs.size(); i++) {
        Vec<N> vec = vecs[i];
        for(size_t j = 0; j < res.size(); j++) {
            vec = vec - res[j] * (g(pt, vecs[i], res[j]) / g(pt, res[j], res[j]) );
        }
        res.push_back(vec);
    }
    if (norm) {
        for(auto& v : res) {
            v = v.normalized();
        }
    }
    return res;
}


template class MetricTensor<1>;
template class MetricTensor<2>;
template class MetricTensor<3>;


template<std::size_t N>
std::vector<Point<N> > geodesic(MetricTensor<N> &g, Point<N> start, Vec<N> vel, std::size_t steps_num, double dt)
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


template std::vector<Point<2> > geodesic(MetricTensor<2> &g, Point<2> start, Vec<2> vel, std::size_t steps_num, double dt);
template std::vector<Point<3> > geodesic(MetricTensor<3> &g, Point<3> start, Vec<3> vel, std::size_t steps_num, double dt);


template<std::size_t N, std::size_t M>
InducedMetricTensor<N, M>::InducedMetricTensor(void(*f)(const Point<N>&, Point<M>&), double pres) :
    f(f),
    pres(pres),
    inv_pres2(0.25 / (pres * pres)),
    inv_pres3(0.125 / (pres * pres * pres))
{
    // This is ALMOST linear version
    this->g = [&, this](Point<N>& p, const Vec<N>& v1, const Vec<N>& v2)
    {
        double alpha1 = v1.norm();
        double alpha2 = v2.norm();

        if(alpha1 * alpha2 == 0) return 0.0; // Almost impossible, but anyway

        this->f(p + v1.normalized() * this->pres, this->dv1p);
        this->f(p - v1.normalized() * this->pres, this->dv1m);
        this->f(p + v2.normalized() * this->pres, this->dv2p);
        this->f(p - v2.normalized() * this->pres, this->dv2m);

        return 0.25 * alpha1 * alpha2 * ((dv1p - dv1m) * (dv2p - dv2m)) / (this->pres * this->pres);
    };
}


template<std::size_t N, std::size_t M>
double InducedMetricTensor<N, M>::dk_gij(Point<N> &pt, size_t k, size_t i, size_t j)
{
    //partialDer<N>(pt, k, [this, i, j](Point<N>& pp){return this->getCoord(pp, i, j);} );
    double pt_k = pt[k];
    double pt_i = pt[i];
    double pt_j = pt[j];

    pt[i] += pres; pt[k] += pres;
    f(pt, vec_pxi_pxk);
    pt[i] = pt_i; pt[k] = pt_k;

    pt[i] -= pres; pt[k] += pres;
    f(pt, vec_mxi_pxk);
    pt[i] = pt_i; pt[k] = pt_k;

    pt[j] += pres; pt[k] += pres;
    f(pt, vec_pxj_pxk);
    pt[j] = pt_j; pt[k] = pt_k;

    pt[j] -= pres; pt[k] += pres;
    f(pt, vec_mxj_pxk);
    pt[j] = pt_j; pt[k] = pt_k;

    pt[i] += pres; pt[k] -= pres;
    f(pt, vec_pxi_mxk);
    pt[i] = pt_i; pt[k] = pt_k;

    pt[i] -= pres; pt[k] -= pres;
    f(pt, vec_mxi_mxk);
    pt[i] = pt_i; pt[k] = pt_k;

    pt[j] += pres; pt[k] -= pres;
    f(pt, vec_pxj_mxk);
    pt[j] = pt_j; pt[k] = pt_k;

    pt[j] -= pres; pt[k] -= pres;
    f(pt, vec_mxj_mxk);
    pt[j] = pt_j; pt[k] = pt_k;


    vec_pxi_pxk.subtract(vec_mxi_pxk);
    vec_pxj_pxk.subtract(vec_mxj_pxk);
    vec_pxi_mxk.subtract(vec_mxi_mxk);
    vec_pxj_mxk.subtract(vec_mxj_mxk);
    double subres1 = vec_pxi_pxk * vec_pxj_pxk;
    double subres2 = vec_pxi_mxk * vec_pxj_mxk;
    /*
    double subres1 = (vec_pxi_pxk - vec_mxi_pxk) * (vec_pxj_pxk - vec_mxj_pxk);
    double subres2 = (vec_pxi_mxk - vec_mxi_mxk) * (vec_pxj_mxk - vec_mxj_mxk);
    */



    return (subres1 - subres2) * inv_pres3;
}

template<std::size_t N, std::size_t M>
double InducedMetricTensor<N, M>::krist(std::size_t l, std::size_t j, std::size_t k, Point<N> &p)
{
    double res = 0;

    Matrix2D<N> Mg;
    for(size_t i = 0; i < N; i++) {
        for(size_t j = 0; j < N; j++) {
            Mg[i][j] = getCoord(p, i, j);
        }
    }
    //Matrix2D<N> MgInv = MetricTensor<N>::getMatrix(p).inverse();
    Matrix2D<N> MgInv = Mg.inverse();

    for(std::size_t r = 0; r < N; r++) {
        double dk_grj = dk_gij(p, k, r, j);
        double dj_grk = dk_gij(p, j, k, r);
        double dr_gjk = dk_gij(p, r, j, k);

        res += MgInv[l][r] * ( dk_grj + dj_grk - dr_gjk );
    }
    return res * 0.5;
}

template<std::size_t N, std::size_t M>
double InducedMetricTensor<N, M>::getCoord(Point<N> &pt, std::size_t i, std::size_t j)
{
    double pt_i = pt[i];
    double pt_j = pt[j];

    Vec<M> vec_pxi, vec_mxi, vec_pxj, vec_mxj;

    pt[i] += pres;
    f(pt, vec_pxi);
    pt[i] = pt_i;

    pt[i] -= pres;
    f(pt, vec_mxi);
    pt[i] = pt_i;

    pt[j] += pres;
    f(pt, vec_pxj);
    pt[j] = pt_j;

    pt[j] -= pres;
    f(pt, vec_mxj);
    pt[j] = pt_j;

    vec_pxi.subtract(vec_mxi);
    vec_pxj.subtract(vec_mxj);
    double subres = vec_pxi * vec_pxj;
    //double subres = (vec_pxi - vec_mxi) * (vec_pxj - vec_mxj);

//    double subres = vec_pxi * vec_pxj -
//                    vec_pxi * vec_mxj -
//                    vec_mxi * vec_pxj +
//                    vec_mxi * vec_mxj;

    return subres * inv_pres2;
}

template class InducedMetricTensor<1,2>;
template class InducedMetricTensor<2,3>;
template class InducedMetricTensor<3,4>;
template class InducedMetricTensor<2,2>;
