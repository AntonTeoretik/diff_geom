#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <iostream>
#include <optional>
#include <algebra.h>

#include <vector>
#include <map>

#include "metric.h"


template <std::size_t N>
using Chart = std::function<bool(Point<N>)>;

template <std::size_t N>
using structMap = std::function<std::optional<Point<N>>(Point<N>)>;

using index = std::size_t;

template <std::size_t N>
struct genPoint {
    index i;
    Point<N> p;
};

template <std::size_t N>
struct genPairOfPoints {
    index i;
    Point<N> p1;
    Point<N> p2;
};

template <class T>
using typedGraph = std::map<std::pair<index, index>, T>;

/**
 * @brief Manifold = Charts + structure maps.
 *
 * The user must ensure that they are diffeomorphism on appropriate domains.
 * @
 */
template <std::size_t N>
class Manifold
{
protected:
    std::vector<Chart<N>> atlas;
    index atlas_size;
    typedGraph<structMap<N>> structureMaps;

public:
    Manifold(const std::vector<Chart<N>>& atlas, const typedGraph<structMap<N>>& structureMaps);

    std::optional<Point<N> > changePointIndex(genPoint<N> pt, index newIndex) const;
    std::optional<Vec<N>> changeVectorIndex(Vec<N> v, genPoint<N> pt, index newIndex, double e=eps) const;
};

template <std::size_t N>
class RiemannianManifold : public Manifold<N>
{
public:
    std::vector<MetricTensor<N>> metric;
    Point<N> doOneStep(Point<N> prev, Point<N> now, index i) const;

public:
    RiemannianManifold(const std::vector<Chart<N>>& atlas,
                       const typedGraph<structMap<N>>& structureMaps,
                       const std::vector<MetricTensor<N>>& metric);

    //pt must be not in the boundary of domain
    std::vector<genPoint<N>> geodesic(genPoint<N> pt, Vec<N> dir, size_t num_of_pts, double step=time_step) const;
};


#endif // MANIFOLD_H
