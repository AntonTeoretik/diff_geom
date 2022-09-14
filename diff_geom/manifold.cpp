#include "manifold.h"


template<std::size_t N>
Manifold<N>::Manifold(const std::vector<Chart<N>>& atlas, const typedGraph<structMap<N>>& structureMaps) :
    atlas(atlas),
    structureMaps(structureMaps)
{

}

template<std::size_t N>
std::optional<genPoint<N> > Manifold<N>::changeIndex(genPoint<N> pt, index newIndex)
{

}

template<std::size_t N>
RiemannianManifold<N>::RiemannianManifold(const std::vector<Chart<N> > &atlas,
                                          const std::vector<std::vector<structMap<N> > > &structureMaps,
                                          const std::vector<MetricTensor<N> > &metric) :
    Manifold<N>(atlas, structureMaps),
    metric(metric)
{

}
