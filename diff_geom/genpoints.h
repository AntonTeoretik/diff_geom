#ifndef GENPOINTS_H
#define GENPOINTS_H

#include "algebra.h"


using chart_index = std::size_t;

template <std::size_t N>
struct genPoint {
    chart_index i;
    Point<N> p;
};

template <std::size_t N, class T>
class genFunction
{
public:
    virtual T operator () (const genPoint<N>&) const = 0;
};

template <std::size_t N>
class diffeomorphism
{
public:
    virtual Point<N> operator () (const Point<N>&) const = 0;
};


#endif // GENPOINTS_H
