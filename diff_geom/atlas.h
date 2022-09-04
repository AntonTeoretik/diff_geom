#ifndef ATLAS_H
#define ATLAS_H

#include <iostream>
#include <optional>
#include <algebra.h>

#include <vector>
#include <map>



template <std::size_t N>
using Chart = std::function<bool(Point<N>)>;

template <std::size_t N>
using structMap = std::function<std::optional<Point<N>>(Point<N>)>;

template <std::size_t N>
using genPoint = std::vector< std::optional<Point<N>> >;

using index = std::size_t;

/**
 * @brief Manifold = Charts + structure maps.
 *
 * The user must ensure that they are diffeomorphism on appropriate domains.
 * @
 */
template <std::size_t N>
class Manifold
{
    std::vector<Chart<N>> atlas;
    std::map<index, std::pair<index, structMap<N>> > structureMaps;
public:
    Manifold(std::vector<Chart<N>> atlas, std::vector<std::vector<structMap<N>>> structureMaps);
};

#endif // ATLAS_H
