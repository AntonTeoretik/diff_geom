#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <iostream>
#include <optional>


#include <vector>
#include <memory>
#include <map>

#include "algebra.h"
#include "metric.h"
#include "genpoints.h"
#include "weights.h"


template <std::size_t N>
class Chart
{
public:
    virtual bool operator ()(const Point<N>&) const = 0;
};

template <class T>
using typedGraph = std::map<std::pair<chart_index, chart_index>, T>;

/**
 * @brief Manifold = Charts + structure maps.
 *
 * The user must ensure that they are diffeomorphism on appropriate domains.
 * @
 */
template <std::size_t N>
class AbstractManifold
{
protected:

    chart_index atlas_size;
public:
    virtual bool changePointIndex(Point<N>& pt, chart_index oldIndex, chart_index newIndex) const = 0;
    virtual bool isPoint(const Point<N>& pt, chart_index i) const = 0;
};

template <std::size_t N>
class AbstractRiemannianManifold : public AbstractManifold<N>
{
public:
    virtual const MetricTensor<N>& getMetric(chart_index i) const = 0;

    Point<N> doOneStep(Point<N> prev, Point<N> now, chart_index i) const;
    void doOneStepWithChange(Point<N>& prev, Point<N>& now, chart_index& i) const;

    std::vector<genPoint<N>> geodesic(const genPoint<N>& pt, Vec<N> dir, size_t num_of_pts, size_t dist = 1, double step=time_step) const;
    std::vector<Vec<N>> orthogonalize(const genPoint<N>& pt, const std::vector<Vec<N> > &vecs, bool normalize = true) const;

    template <class T>
    T integrateAlongPath(genPoint<N> start,
                     Vec<N> dir,
                     size_t num_of_pts,
                     const genFunction<N,T>& func,
                     const weightFunction& weight,
                     double step=time_step) const;
};

template <class T, std::size_t N>
T integrateAlongPath(const std::vector<genPoint<N>>& points,
                     const std::function<T(genPoint<N>)>& func,
                     const std::function<double(double)>& weight=[](double){return 1.0;},
                     double step=time_step);


#endif // MANIFOLD_H
