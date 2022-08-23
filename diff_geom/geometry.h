#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <functional>
#include <array>
#include <constants.h>

#include <algebra.h>



template<std::size_t R, std::size_t N>
class AbstractTensor
{
public:
    virtual double operator() (Point<N>, std::array<Vec<N>, R>) = 0;
    // In each implementation the user must be sure that it is a polylinear form in each point!
};

template<std::size_t N>
class MetricTensor : public AbstractTensor<2,N>
{
    std::function<double(Point<N>, Vec<N>, Vec<N>)> g; // Must be symmetric, non degenerate
public:
    MetricTensor<N>();
    MetricTensor<N>(std::function<double(Point<N>, Vec<N>, Vec<N>)>);


    virtual double operator() (Point<N>, std::array<Vec<N>, 2>);
    virtual double operator() (Point<N>, Vec<N>, Vec<N>);

    double krist(std::size_t l, std::size_t j , std::size_t k, Point<N> p);
};



#endif // GEOMETRY_H
