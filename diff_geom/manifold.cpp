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
std::optional<Vec<N> > Manifold<N>::changeVectorIndex(Vec<N> v, genPoint<N> pt, index newIndex, double eps)
{
    if(v == Vec<N>{0.0,0.0}) return {{0.0,0.0}};
    auto norm = v.norm();
    auto nv = v.nomalized();
    auto p = pt.p;

    if( auto pp_next = changePointIndex({pt.i, p + nv * eps}, newIndex),
        pp_prev = changePointIndex({pt.i, p - nv * eps}, newIndex);
        pp_next.has_value() and pp_prev.has_value()) {

        return { (pp_next.value().p - pp_prev.value().p) * (norm*0.5/eps) };
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
