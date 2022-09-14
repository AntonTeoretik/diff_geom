#include "manifold.h"


template<std::size_t N>
Manifold<N>::Manifold(const std::vector<Chart<N>>& atlas, const typedGraph<structMap<N>>& structureMaps) :
    atlas(atlas),
    atlas_size(atlas.size()),
    structureMaps(structureMaps)
{

}

template<std::size_t N>
std::optional<genPoint<N> > Manifold<N>::changePointIndex(genPoint<N> pt, index newIndex)
{
    if (newIndex >= atlas_size) return {};
    if (auto it = structureMaps.find({pt.i, newIndex}); it != structureMaps.end())
    {
        if(auto p = (it->second)(pt.p); p) return { {newIndex, p.value()} };
        return {};
    }
    return {};
}

template<std::size_t N>
RiemannianManifold<N>::RiemannianManifold(const std::vector<Chart<N> > &atlas,
                                          const std::vector<std::vector<structMap<N> > > &structureMaps,
                                          const std::vector<MetricTensor<N> > &metric) :
    Manifold<N>(atlas, structureMaps),
    metric(metric)
{

}


template class Manifold<2>;
template class Manifold<3>;
