#ifndef METRIC_H
#define METRIC_H

#include <functional>
#include <array>
#include <constants.h>

#include <vector>
#include <algebra.h>



template<std::size_t R, std::size_t N>
class AbstractTensor
{
public:
    virtual double operator() (Point<N>&, std::array<Vec<N>, R>&) const = 0;
    // In each implementation the user must be sure that it is a polylinear form in each point!
};

template<std::size_t N>
class MetricTensor : public AbstractTensor<2,N>
{
protected:
    std::function<double(Point<N>&, const Vec<N>&, const Vec<N>&)> g; // Must be symmetric, non degenerate
public:
    MetricTensor<N>();
    MetricTensor<N>(std::function<double(Point<N>&, const Vec<N>&, const Vec<N>&)>);

    virtual double operator() (Point<N>&, std::array<Vec<N>, 2>&) const;
    virtual double operator() (Point<N>&, Vec<N>&, Vec<N>&) const;

    virtual Matrix2D<N> getMatrix(Point<N>& P) const;
    virtual double getCoord(Point<N>& P, std::size_t i, std::size_t j) const;

    double krist(std::size_t l, std::size_t j , std::size_t k, Point<N>& p) const; // Returns Г^l_jk
    Matrix2D<N> kristMatrix(std::size_t l, Point<N>& p) const; // Returns Г^l_{jk}

    std::vector<Vec<N>> orthogonalize(Point<N>& pt, const std::vector<Vec<N> > &vecs, bool normalize=true) const;
};

template <std::size_t N, std::size_t M>
class InducedMetricTensor : public MetricTensor<N>
{
protected:
    std::function<Point<M>(const Point<N>&)> f;
    double dk_gij(Point<N>& pt, size_t k, size_t r, size_t j) const;
    double pres;

public:
    InducedMetricTensor<N, M>(std::function<Point<M>(const Point<N>&)>, double pres=eps);
    double krist(std::size_t l, std::size_t j , std::size_t k, Point<N>& p) const; // Returns Г^l_jk

    Point<M> apply_generator(const Point<N>&) const;
};


template<std::size_t N>
inline double partialDer(const Point<N>& A, Vec<N> dir, const std::function<double(const Point<N>&)>& f, double step = eps)
{
    return (f(A + dir * step) - f(A + dir * (-step))) / (2 * step) ;
}

//For basis vector
template<std::size_t N>
inline double partialDer(Point<N>& A, size_t i, const std::function<double(Point<N>&)>& f, double step = eps)
{
    double savedValue = A[i];
    A[i] += step;
    double f_next = f(A);
    A[i] = savedValue - step;
    double f_prev = f(A);
    A[i] = savedValue;

    return (f_next - f_prev) / (2 * step);
}

// Geodesic 1

template<std::size_t N>
std::vector<Point<N>> geodesic(const MetricTensor<N> &g, Point<N> start, Vec<N> vel, std::size_t steps_num, double dt = time_step);


#endif // METRIC_H
