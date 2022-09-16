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
    //std::cout << "changePointIndex: " << pt.i << " " << pt.p.to_str() << std::endl;
    if (newIndex >= atlas_size) return {};



    if (auto it = structureMaps.find({pt.i, newIndex}); it != structureMaps.end())
    {
        if(auto p = (it->second)(pt.p); p.has_value()) return { p.value() };
        return {};
    }
    return {};
}

template<std::size_t N>
std::optional<Vec<N> > Manifold<N>::changeVectorIndex(Vec<N> v, genPoint<N> pt, index newIndex, double eps) const
{
    if(v == Vec<N>::zero()) return {Vec<N>::zero()};
    auto norm = v.norm();
    auto nv = v.normalized();
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

    Vec<N> curv;
    for (size_t l = 0; l < N; l++) {
        curv[l] = (now - prev) * (metric[i]->kristMatrix(l, now) * (now - prev));
    }


    auto next = now * 2 - prev - curv;

    return next;

}

template<std::size_t N>
RiemannianManifold<N>::RiemannianManifold(const std::vector<Chart<N> > &atlas,
                                          const typedGraph<structMap<N> > &structureMaps,
                                          const std::vector<std::shared_ptr<MetricTensor<N> > > &metric) :
    Manifold<N>(atlas, structureMaps),
    metric(metric)
{

}

template<std::size_t N>
std::vector<genPoint<N> > RiemannianManifold<N>::geodesic(genPoint<N> pt, Vec<N> dir, size_t num_of_pts, size_t dist, double step) const
{
    Point<N> prev = pt.p;
    Point<N> now = prev + (dir * step); // Possible trouble here
    Point<N> next;

    size_t count = 1;

    std::vector<genPoint<N> > res = {pt};
    if (dist == 1) {
        res.push_back({pt.i, now});
        count = 0;
    }

    index cur_index = pt.i;

    for (size_t i = 0; i < num_of_pts; i++) {
        count++;

        next = doOneStep(prev, now, cur_index);

        // next is in current domain
        if(this->atlas[cur_index](next)) {
            prev = now;
            now = next;

            if (count == dist) {
                count = 0;
                res.push_back({cur_index, now});
            }
            //std::cout << cur_index << ": " << now.to_str() << std::endl;
            continue;
        }
        // next point is outside -> need to change domain.
        else {
            bool change_failed = true;
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

                        if (count == dist) {
                            count = 0;

                            res.push_back({cur_index, now});
                        }

                        change_failed = false;
                        break;
                    }
                }
            }
            if(change_failed) {
                std::cout << "Failed to change the index: cur_index = " << cur_index << ": " << now.to_str() << ", " << next.to_str()  << std::endl;
            }
        }
        // This situation must be impossible!
    }
    return res;
}

template<std::size_t N>
std::vector<Vec<N> > RiemannianManifold<N>::orthogonalize(genPoint<N> pt, const std::vector<Vec<N> > &vecs, bool normalize) const
{
    return this->metric[pt.i]->orthogonalize(pt.p, vecs, normalize);
}


template class RiemannianManifold<1>;
template class RiemannianManifold<2>;
template class RiemannianManifold<3>;


template<class T, std::size_t N>
T integrateAlongPath(const std::vector<genPoint<N> > &points,
                     const std::function<T (genPoint<N>)> func,
                     const std::function<double (double)> weight,
                     double step)
{
    double current_time = 0.0;
    T res = {};
    //first value
    if(points.size() == 0) {
        return res;
    }
    res = func(points[0]) * (step * weight(current_time));

    for(size_t i = 1; i <= points.size(); i++) {
        current_time += step;
        res = res + func(points[0]) * (step * weight(current_time));
    }
    return res;
}

template
double integrateAlongPath(const std::vector<genPoint<3> > &,
                          const std::function<double (genPoint<3>)>,
                          const std::function<double (double)>,
                          double);

template
double integrateAlongPath(const std::vector<genPoint<2> > &,
                          const std::function<double (genPoint<2>)>,
                          const std::function<double (double)>,
                          double);

template
Vec<3> integrateAlongPath(const std::vector<genPoint<3> > &,
                          const std::function<Vec<3> (genPoint<3>)>,
                          const std::function<double (double)>,
                          double);
