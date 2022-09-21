#ifndef METRIC_H
#define METRIC_H

#include <functional>
#include <array>
#include <vector>

#include "constants.h"
#include "algebra.h"



template<std::size_t R, std::size_t N>
class AbstractTensor
{
public:
    virtual double operator() (const Point<N>&, const std::array<Vec<N>, R>&) const = 0;
    // In each implementation the user must be sure that it is a polylinear form in each point!
};





template<std::size_t N>
class MetricTensor : public AbstractTensor<2,N>
{
protected:
public:
    MetricTensor<N>();

    virtual double operator() (const Point<N>&, const std::array<Vec<N>, 2>&) const;
    virtual double operator() (const Point<N>&, const Vec<N>&, const Vec<N>&) const = 0;

    virtual Matrix2D<N> getMatrix(const Point<N>& P) const;
    virtual double getCoord(const Point<N>& P, std::size_t i, std::size_t j) const;

    virtual double krist(std::size_t l, std::size_t j , std::size_t k, const Point<N>& p) const; // Returns Г^l_jk
    Matrix2D<N> kristMatrix(std::size_t l, const Point<N>& p) const; // Returns Г^l_{jk}

    std::vector<Vec<N>> orthogonalize(const Point<N>& pt, const  std::vector<Vec<N> > &vecs, bool normalize=true) const ;
};



template <std::size_t N, std::size_t M>
class InducedMetricTensor : public MetricTensor<N>
{
protected:    
    double dk_gij(const Point<N>& pt, size_t k, size_t r, size_t j) const;
    double pres, inv_pres2, inv_pres3;

    Vec<M> dv1p, dv1m, dv2p, dv2m;
    Vec<M> vec_pxi, vec_mxi, vec_pxj, vec_mxj;

public:
    InducedMetricTensor<N, M>(double pres=eps);

    virtual Point<M> gen_func(const Point<N>&) const = 0;

    virtual double operator() (const Point<N>&, const Vec<N>&, const Vec<N>&) const;

    virtual double krist(std::size_t l, std::size_t j , std::size_t k, const Point<N>& p) const; // Returns Г^l_jk
    virtual double getCoord(const Point<N>& P, std::size_t i, std::size_t j) const;
};


//For basis vector
template<std::size_t N>
inline double partialDer(const Point<N>& A_, size_t i, const std::function<double(const Point<N>&)>& f, double step = eps)
{
    Point<N> A = A_;
    double savedValue = A[i];
    A[i] += step;
    double f_next = f(A);
    A[i] = savedValue - step;
    double f_prev = f(A);
    A[i] = savedValue;

    return (f_next - f_prev) / (2 * step);
}

//template<std::size_t N>
//inline double partialDer(const Point<N>& A, const Vec<N>& dir,  std::function<double(const Point<N>&)>& f, double step = eps)
//{
//    return (f(A + dir * step) - f(A + dir * (-step))) / (2 * step) ;
//}





#endif // METRIC_H
