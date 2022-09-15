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
    virtual double operator() (Point<N>, std::array<Vec<N>, R>) const = 0;
    // In each implementation the user must be sure that it is a polylinear form in each point!
};

template<std::size_t N>
class MetricTensor : public AbstractTensor<2,N>
{
protected:
    std::function<double(Point<N>, Vec<N>, Vec<N>)> g; // Must be symmetric, non degenerate
public:
    MetricTensor<N>();
    MetricTensor<N>(std::function<double(Point<N>, Vec<N>, Vec<N>)>);

    virtual double operator() (Point<N>, std::array<Vec<N>, 2>) const;
    virtual double operator() (Point<N>, Vec<N>, Vec<N>) const;

    virtual Matrix2D<N> getMatrix(Point<N> P) const;
    virtual double getCoord(Point<N> P, std::size_t i, std::size_t j) const;

    double krist(std::size_t l, std::size_t j , std::size_t k, Point<N> p) const; // Returns Г^l_jk
    Matrix2D<N> kristMatrix(std::size_t l, Point<N> p) const; // Returns Г^l_{jk}
};

template <std::size_t N, std::size_t M>
class InducedMetricTensor : public MetricTensor<N>
{
protected:
    std::function<Point<M>(Point<N>)> f;
    double pres;

public:
    InducedMetricTensor<N, M>(std::function<Point<M>(Point<N>)>, double pres=eps);

    Point<M> apply_generator(Point<N>) const;
};

template<std::size_t N>
inline double partialDer(Point<N> A, Vec<N> dir, std::function<double(Point<N>)> f, double step = eps)
{
    return (f(A + dir * step) - f(A + dir * (-step))) / (2 * step) ;
}


// Geodesic 1

template<std::size_t N>
std::vector<Point<N>> geodesic(const MetricTensor<N> &g, Point<N> start, Vec<N> vel, std::size_t steps_num, double dt = time_step);


#endif // METRIC_H
