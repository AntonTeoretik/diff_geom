#include "metric.h"



template<std::size_t N>
CUDA_F MetricTensor<N>::MetricTensor()
{
    for(std::size_t i = 0; i < N; i++) {
        for(std::size_t j = 0; j < N; j++) {
            _basis[i][j] = (i == j ? 1.0 : 0.0);
        }
    }
}

template<std::size_t N>
CUDA_F double MetricTensor<N>::call_array(const Point<N>& p, const Vec<N>* vecs) const
{
    return (*this)(p, vecs[0], vecs[1]);
}


template<std::size_t N>
CUDA_F Matrix2D<N> MetricTensor<N>::getMatrix(const Point<N>& P) const
{
    Matrix2D<N> m;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            m.get(i, j) = getCoord(P, i, j);
        }
    }
    return m;
}

template<std::size_t N>
CUDA_F double MetricTensor<N>::getCoord(const Point<N>& P, std::size_t i, std::size_t j) const
{
    // printf("GetCoord %lu %lu %lu\n", N, i, j);
    // return 0.0;
    return (*this)(P, _basis[i], _basis[j]);
}

// template<std::size_t N>
// CUDA_F double MetricTensor<N>::krist(std::size_t l, std::size_t j, std::size_t k, const Point<N>& p) const
// {
//     double res = 0;
//     return 0;

//     // Matrix2D<N> MgInv = getMatrix(p).inverse();

//     // for(std::size_t r = 0; r < N; r++) {
//     //     double dk_grj = partialDer<N>(p, k, [this, r, j] __device__ (const Point<N>& pp){return this->getCoord(pp, r, j);} );
//     //     double dj_grk = partialDer<N>(p, j, [this, r, k] __device__ (const Point<N>& pp){return this->getCoord(pp, r, k);} );
//     //     double dr_gjk = partialDer<N>(p, r, [this, k, j] __device__ (const Point<N>& pp){return this->getCoord(pp, j, k);} );

//     //     res += MgInv.get(l, r) * ( dk_grj + dj_grk - dr_gjk );
//     // }
//     // return res * 0.5;
// }

template<std::size_t N>
CUDA_F Matrix2D<N> MetricTensor<N>::kristMatrix(std::size_t l, const Point<N>& p) const
{
    Matrix2D<N> res;
    // return res;
    for (size_t i = 1; i < N; i++) {
        for (size_t j = 0; j < i; j++) {
            double kr = krist(l, i, j, p);
            res.get(i, j) = kr;
            res.get(j, i) = kr;
        }
    }
    for (size_t i = 0; i < N; i++) {
        res.get(i, i) = this->krist(l, i, i, p);
    }

    //return Matrix2D<N>([this, l, &p](auto i, auto j){return this->krist(l, i, j, p);});
    return res;
}


template<std::size_t N>
CUDA_F void MetricTensor<N>::orthogonalize(
    const Point<N>& pt, const Span<Vec<N>>& vecs, bool norm, 
    Span<Vec<N>>& res, int& n_out) const
{
    n_out = 0;
    // std::vector<Vec<N> > res = {};
    for (size_t i = 0; i < vecs.size; i++) {
        Vec<N> vec = vecs.data[i];
        for(size_t j = 0; j < n_out; j++) {
            vec = vec - res[j] * ((*this)(pt, vecs[i], res[j]) / (*this)(pt, res[j], res[j]) );
        }
        res.data[n_out] = vec;
        ++n_out;
    }
    if (norm) {
        for(int i = 0; i < n_out; ++i) {
            //std::cout << "Scale: " << std::sqrt((*this)(pt, v, v)) << std::endl;
            //std::cout << "v: " << v.to_str() << std::endl;
            auto v = res.data[i];
            res.data[i].scale(1.0 / std::sqrt((*this)(pt, v, v)));
        }
    }
}


template class MetricTensor<1>;
template class MetricTensor<2>;
template class MetricTensor<3>;


template<std::size_t N, std::size_t M>
CUDA_F InducedMetricTensor<N, M>::InducedMetricTensor(double pres) :
    MetricTensor<N>(),
    pres(pres),
    inv_pres2(0.25 / (pres * pres)),
    inv_pres3(0.125 / (pres * pres * pres))
{
    //std::cout << pres << " " << inv_pres2 << " " << inv_pres3 << std::endl;
}


template<std::size_t N, std::size_t M>
CUDA_F double InducedMetricTensor<N, M>::operator()(const Point<N> & p, const Vec<N> & v1_, const Vec<N> & v2_) const
{
    Vec<N> v1 = v1_;
    Vec<N> v2 = v2_;
    double alpha1 = v1.norm();
    double alpha2 = v2.norm();
    v1.normalize();
    v2.normalize();
    v1.scale(pres);
    v2.scale(pres);

    Vec<M> dv1p, dv1m, dv2p, dv2m;

    if(alpha1 * alpha2 == 0) return 0.0; // Almost impossible, but anyway

    dv1p = gen_func(p + v1);
    dv1m = gen_func(p - v1);
    dv2p = gen_func(p + v2);
    dv2m = gen_func(p - v2);
    dv1p.subtract(dv1m);
    dv2p.subtract(dv2m);

    return alpha1 * alpha2 * (dv1p * dv2p) * 0.25 / pres / pres;
};

template<std::size_t N, std::size_t M>
CUDA_F double InducedMetricTensor<N, M>::dk_gij(const Point<N> &p, size_t k, size_t i, size_t j) const
{
    //partialDer<N>(pt, k, [this, i, j](Point<N>& pp){return this->getCoord(pp, i, j);} );

    Point<N> pt = p;

    double pt_k = pt[k];
    double pt_i = pt[i];
    double pt_j = pt[j];

    Vec<M> vec_pxi, vec_mxi, vec_pxj, vec_mxj;

    pt[i] += pres; pt[k] += pres;
    vec_pxi = gen_func(pt);
    pt[i] = pt_i; pt[k] = pt_k;

    pt[i] -= pres; pt[k] += pres;
    vec_mxi = gen_func(pt);
    pt[i] = pt_i; pt[k] = pt_k;

    pt[j] += pres; pt[k] += pres;
    vec_pxj = gen_func(pt);
    pt[j] = pt_j; pt[k] = pt_k;

    pt[j] -= pres; pt[k] += pres;
    vec_mxj = gen_func(pt);
    pt[j] = pt_j; pt[k] = pt_k;

    vec_pxi.subtract(vec_mxi);
    vec_pxj.subtract(vec_mxj);

    double subres1 = vec_pxi * vec_pxj;


    pt[i] += pres; pt[k] -= pres;
    vec_pxi = gen_func(pt);
    pt[i] = pt_i; pt[k] = pt_k;

    pt[i] -= pres; pt[k] -= pres;
    vec_mxi = gen_func(pt);
    pt[i] = pt_i; pt[k] = pt_k;

    pt[j] += pres; pt[k] -= pres;
    vec_pxj = gen_func(pt);
    pt[j] = pt_j; pt[k] = pt_k;

    pt[j] -= pres; pt[k] -= pres;
    vec_mxj = gen_func(pt);
    pt[j] = pt_j; pt[k] = pt_k;

    vec_pxi.subtract(vec_mxi);
    vec_pxj.subtract(vec_mxj);

    double subres2 = vec_pxi * vec_pxj;

    return (subres1 - subres2) * 0.125 / pres / pres / pres;
}

template<std::size_t N, std::size_t M>
CUDA_F double InducedMetricTensor<N, M>::krist(std::size_t l, std::size_t j, std::size_t k, const Point<N> &p) const
{
    double res = 0;
    // double x = getCoord(p, 0, 0);;
    // return res;

    Matrix2D<N> Mg;
    for(size_t i = 0; i < N; i++) {
        for(size_t jj = 0; jj < N; jj++) {
            Mg.get(i, jj) = getCoord(p, i, jj);
        }
    }
    //Matrix2D<N> MgInv = MetricTensor<N>::getMatrix(p).inverse();
    Matrix2D<N> MgInv = Mg.inverse();

    for(std::size_t r = 0; r < N; r++) {
        double dk_grj = dk_gij(p, k, r, j);
        double dj_grk = dk_gij(p, j, k, r);
        double dr_gjk = dk_gij(p, r, j, k);

        res += MgInv.get(l, r) * ( dk_grj + dj_grk - dr_gjk );
    }
    return res * 0.5;
}

template<std::size_t N, std::size_t M>
CUDA_F double InducedMetricTensor<N, M>::getCoord(const Point<N> &p, std::size_t i, std::size_t j) const
{
    // printf("GetCoord\n");
    Point<N> pt = p;
    double pt_i = pt[i];
    double pt_j = pt[j];

    Vec<M> vec_pxi, vec_mxi, vec_pxj, vec_mxj;

    pt[i] += pres;
    vec_pxi = gen_func(pt);
    pt[i] = pt_i;

    pt[i] -= pres;
    vec_mxi = gen_func(pt);
    pt[i] = pt_i;

    pt[j] += pres;
    vec_pxj = gen_func(pt);
    pt[j] = pt_j;

    pt[j] -= pres;
    vec_mxj = gen_func(pt);
    pt[j] = pt_j;

    vec_pxi.subtract(vec_mxi);
    vec_pxj.subtract(vec_mxj);
    double subres = vec_pxi * vec_pxj;
    return 0.25 * subres / pres / pres;
}

template class InducedMetricTensor<1,2>;
template class InducedMetricTensor<2,3>;
template class InducedMetricTensor<3,4>;
template class InducedMetricTensor<2,2>;
