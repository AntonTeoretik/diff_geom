#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <iostream>
#include <optional>
#include <algebra.h>


#include <vector>
#include <memory>
#include <map>

#include "metric.h"


template <std::size_t N>
using Chart = std::function<bool(Point<N>)>;

template <std::size_t N>
using structMap = std::function<std::optional<Point<N>>(Point<N>)>;

using chart_index = std::size_t;

template <std::size_t N>
struct genPoint {
    chart_index i;
    Point<N> p;
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
class Manifold
{
public:
    std::vector<Chart<N>> atlas;
    chart_index atlas_size;
    typedGraph<structMap<N>> structureMaps;

public:
    Manifold(const std::vector<Chart<N>>& atlas, const typedGraph<structMap<N>>& structureMaps);

    std::optional<Point<N> > changePointIndex(genPoint<N> pt, chart_index newIndex) const;
    std::optional<Vec<N>> changeVectorIndex(Vec<N> v, genPoint<N> pt, chart_index newIndex, double e=eps) const;
};

template <std::size_t N>
class RiemannianManifold : public Manifold<N>
{
public:
    std::vector<std::shared_ptr<MetricTensor<N>>> metric; //TODO
    Point<N> doOneStep(Point<N> prev, Point<N> now, chart_index i) const;

public:
    RiemannianManifold(const std::vector<Chart<N>>& atlas,
                       const typedGraph<structMap<N>>& structureMaps,
                       const std::vector<std::shared_ptr<MetricTensor<N>>>& metric);

    //pt must be not in the boundary of domain
    std::vector<genPoint<N>> geodesic(genPoint<N> pt, Vec<N> dir, size_t num_of_pts, size_t dist = 1, double step=time_step) const;
    std::vector<Vec<N>> orthogonalize(genPoint<N> pt, const std::vector<Vec<N> > &vecs, bool normalize = true) const;
};


template <class T, std::size_t N>
T integrateAlongPath(const std::vector<genPoint<N>>& points,
                     const std::function<T(genPoint<N>)> func,
                     const std::function<double(double)> weight=[](double){return 1.0;},
                     double step=time_step);


#endif // MANIFOLD_H
