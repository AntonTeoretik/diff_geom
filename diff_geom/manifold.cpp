#include "manifold.h"


template<std::size_t N>
Manifold<N>::Manifold(const std::vector<Chart<N>>& atlas, const typedGraph<structMap<N>>& structureMaps) :
    atlas(atlas),
    atlas_size(atlas.size()),
    structureMaps(structureMaps)
{

}

template<std::size_t N>
std::optional<Point<N> > Manifold<N>::changePointIndex(genPoint<N> pt, index newIndex) const
{
    if (newIndex >= atlas_size) return {};
    if (auto it = structureMaps.find({pt.i, newIndex}); it != structureMaps.end())
    {
        if(auto p = (it->second)(pt.p); p) return { p.value() };
        return {};
    }
    return {};
}

template<std::size_t N>
std::optional<Vec<N> > Manifold<N>::changeVectorIndex(Vec<N> v, genPoint<N> pt, index newIndex, double eps) const
{
    if(v == Vec<N>::zero()) return {Vec<N>::zero()};
    auto norm = v.norm();
    auto nv = v.nomalized();
    auto p = pt.p;

    if( auto pp_next = changePointIndex({pt.i, p + nv * eps}, newIndex),
        pp_prev = changePointIndex({pt.i, p - nv * eps}, newIndex);
        pp_next.has_value() and pp_prev.has_value()) {

        return { (pp_next.value() - pp_prev.value()) * (norm*0.5/eps) };
    }
    return {};
}

template class Manifold<1>;
template class Manifold<2>;
template class Manifold<3>;


template<std::size_t N>
Point<N> RiemannianManifold<N>::doOneStep(Point<N> prev, Point<N> now, index i) const
{
    //Try to define next point in this domain.
    std::cout << "DoOneStep? " << i << std::endl;
    Vec<N> curv;
    for (size_t l = 0; l < N; l++) {
        auto krMat = metric[i].kristMatrix(l, now);

        std::cout << krMat.to_str() << std::endl;

        curv[l] = (now - prev) * (metric[i].kristMatrix(l, now) * (now - prev));
        std::cout << "DoOneStep? " << i << " " << l << std::endl;
    }


    auto next = now * 2 - prev - curv;

    return next;

}

template<std::size_t N>
RiemannianManifold<N>::RiemannianManifold(const std::vector<Chart<N> > &atlas,
                                          const typedGraph<structMap<N> > &structureMaps,
                                          const std::vector<MetricTensor<N> > &metric) :
    Manifold<N>(atlas, structureMaps),
    metric(metric)
{

}

template<std::size_t N>
std::vector<genPoint<N> > RiemannianManifold<N>::geodesic(genPoint<N> pt, Vec<N> dir, size_t num_of_pts, double step) const
{
    Point<N> prev = pt.p;
    Point<N> now = prev + (dir * step); // Possible trouble here
    Point<N> next;
    std::cout << "Hey" << std::endl;

    std::vector<genPoint<N> > res = {pt, {pt.i, now}};

    index cur_index = pt.i;

    for (size_t i = 0; i < num_of_pts; i++) {

        next = doOneStep(prev, now, cur_index);

        std::cout << "Hey: " << i << std::endl;

        // next is in current domain
        if(this->atlas[cur_index](next)) {
            prev = now;
            now = next;
            res.push_back({cur_index, now});
            continue;
        }
        // next point is outside -> need to change domain.
        for(index new_index = 0; new_index < this->atlas_size; new_index++) {
            if (std::optional<Point<N>> alt_prev_op = this->changePointIndex({cur_index, prev}, new_index),
                                        alt_now_op  = this->changePointIndex({cur_index, now}, new_index);
                alt_prev_op.has_value() and alt_now_op.has_value())
            {
                auto alt_next = doOneStep(alt_prev_op.value(), alt_now_op.value(), new_index);
                // Everything is good
                if(this->atlas[new_index](alt_next)) {
                    prev = alt_now_op.value();
                    now = alt_next;
                    cur_index = new_index;

                    res.push_back({cur_index, now});
                    break;
                }
            }
        }
        // This situation must be impossible!
    }
    return res;
}


template class RiemannianManifold<1>;
template class RiemannianManifold<2>;
template class RiemannianManifold<3>;

