#ifndef METRIC_H
#define METRIC_H

#include <functional>
#include <array>

#include "constants.h"
#include "algebra.h"
#include "util.h"


template<std::size_t R, std::size_t N>
class AbstractTensor
{
public:
    virtual CUDA_F double call_array (const Point<N>&, const Vec<N>*) const = 0;
    // In each implementation the user must be sure that it is a polylinear form in each point!
};

template<std::size_t N>
class MetricTensor : public AbstractTensor<2,N>
{
protected:
    Vec<N> _basis[N];
public:
    CUDA_F MetricTensor<N>();

    virtual CUDA_F double call_array (const Point<N>&, const Vec<N>*) const override;
    virtual CUDA_F double operator() (const Point<N>&, const Vec<N>&, const Vec<N>&) const = 0;

    virtual CUDA_F Matrix2D<N> getMatrix(const Point<N>& P) const;
    virtual CUDA_F double getCoord(const Point<N>& P, std::size_t i, std::size_t j) const;

    virtual CUDA_F double krist(std::size_t l, std::size_t j , std::size_t k, const Point<N>& p) const = 0; // Returns Г^l_jk
    CUDA_F Matrix2D<N> kristMatrix(std::size_t l, const Point<N>& p) const; // Returns Г^l_{jk}

    CUDA_F void orthogonalize(
        const Point<N>& pt, const Span<Vec<N>> &vecs, bool normalize,
        Span<Vec<N>>& vecs_out, int& n_out) const ;
};



template <std::size_t N, std::size_t M>
class InducedMetricTensor : public MetricTensor<N>
{
protected:    
    double CUDA_F dk_gij(const Point<N>& pt, size_t k, size_t r, size_t j) const;
    double pres, inv_pres2, inv_pres3;

public:
    CUDA_F InducedMetricTensor<N, M>(double pres=eps);

    virtual CUDA_F Point<M> gen_func(const Point<N>&) const = 0;

    virtual CUDA_F double operator() (const Point<N>&, const Vec<N>&, const Vec<N>&) const override;

    virtual CUDA_F double krist(std::size_t l, std::size_t j , std::size_t k, const Point<N>& p) const override; // Returns Г^l_jk
    virtual CUDA_F double getCoord(const Point<N>& P, std::size_t i, std::size_t j) const override;
};


//For basis vector
template<std::size_t N, typename Func>
inline CUDA_F double partialDer(const Point<N>& A_, size_t i, Func f, double step = eps)
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

template<std::size_t N, typename Func>
inline CUDA_F double partialDer(const Point<N>& A, const Vec<N>& dir, Func f, double step = eps)
{
    return (f(A + dir * step) - f(A + dir * (-step))) / (2 * step) ;
}




#endif // METRIC_H
