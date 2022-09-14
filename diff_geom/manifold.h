#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <iostream>
#include <optional>
#include <algebra.h>

#include <vector>
#include <map>

#include "geometry.h"


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

    std::optional<genPoint<N>> changePointIndex(genPoint<N> pt, index newIndex);
    std::optional<Vec<N>> changeVectorIndex(Vec<N> v, genPoint<N> pt, index newIndex, double e=eps);
};

template <std::size_t N>
class RiemannianManifold : public Manifold<N>
{
protected:
    std::vector<MetricTensor<N>> metric;
public:
    RiemannianManifold(const std::vector<Chart<N>>& atlas,
                       const std::vector<std::vector<structMap<N>>>& structureMaps,
                       const std::vector<MetricTensor<N>>& metric);

};


#endif // MANIFOLD_H
