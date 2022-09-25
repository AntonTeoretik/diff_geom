#include "manifold.h"



// template<std::size_t N>
// CUDA_F void AbstractRiemannianManifold<N>::geodesic(
//         const genPoint<N> &pt,
//         Vec<N> dir,
//         size_t num_of_pts,
//         size_t dist,
//         double step,
//         Span<genPoint<N>>& res,
//         int& n_out
//         ) const
// {
//     Point<N> prev = pt.p;
//     Point<N> now = prev + (dir * step); // Possible trouble here

//     size_t count = 1;

//     res.data[0] = pt;
//     n_out = 1;
//     if (dist == 1) {
//         res.data[n_out] = genPoint<N>{pt.i, now};
//         count = 0;
//         ++n_out;
//     }

//     chart_index cur_index = pt.i;

//     for (size_t i = 0; i < num_of_pts; i++) {
//         count++;
//         doOneStepWithChange(prev, now, cur_index);

//         if (count == dist) {
//             count = 0;
//             res.data[n_out] = getPoint<N>{cur_index, now};
//             ++n_out;
//         }
//     }
// }


template<std::size_t N>
CUDA_F void AbstractRiemannianManifold<N>::orthogonalize(
    const genPoint<N> &pt, const Span<Vec<N> > &vecs, bool normalize,
    Span<Vec<N>>& out, int& n_out) const
{
    return this->getMetric(pt.i).orthogonalize(pt.p, vecs, normalize, out, n_out);
}


// template<class T, std::size_t N>
// T integrateAlongPath(const std::vector<genPoint<N> > &points,
//                      const std::function<T (genPoint<N>)> &func,
//                      const std::function<double (double)> &weight,
//                      double step)
// {
//     double current_time = 0.0;
//     T res = {};
//     //first value
//     if(points.size() == 0) {
//         return res;
//     }
//     res = func(points[0]) * (step * weight(current_time));

//     for(size_t i = 1; i <= points.size(); i++) {
//         current_time += step;
//         res = res + func(points[i]) * (step * weight(current_time));
//     }
//     return res;
// }

template<size_t N>
template<class T>
CUDA_F T AbstractRiemannianManifold<N>::integrateAlongPath(genPoint<N> start,
                                                  Vec<N> dir, size_t num_of_pts,
                                                  const genFunction<N, T> &func,
                                                  const weightFunction &weight,
                                                  double step) const
{
    //std::cout << "AbstractRiemannianManifold<N>::integrateAlongPath : " << dir.to_str() << std::endl;

    double current_time = 0.0;

    Point<N> prev = start.p;
    Point<N> now = prev + (dir * step); // Possible trouble here

    chart_index cur_index = start.i;

    //NORMALIZE DIR
    double norm2 = getMetric(cur_index)(prev, dir, dir);
    dir.scale(1.0 / sqrt(norm2));

    T res = func(genPoint<N>{cur_index, now}) * weight(0.0) * step;

    for(size_t i = 1; i <= num_of_pts; i++) {
        current_time += step;
        doOneStepWithChange(prev, now, cur_index);
        genPoint<N> pt{cur_index, now};
        res = res + func(pt) * (step * weight(current_time));
    }
    return res;
}

template<std::size_t N>
Point<N> AbstractRiemannianManifold<N>::doOneStep(Point<N> prev, Point<N> now, chart_index i) const
{
    //Try to define next point in this domain.

    Vec<N> curv;
    for (size_t l = 0; l < N; l++) {
        curv[l] = (now - prev) * (getMetric(i).kristMatrix(l, now) * (now - prev));
    }

    auto next = now * 2 - prev - curv;

    return next;
}

template<std::size_t N>
CUDA_F void AbstractRiemannianManifold<N>::doOneStepWithChange(Point<N> &prev, Point<N> &now, chart_index cur_index) const
{
    //std::cout << prev.to_str() << "; " << now.to_str() << "; " << cur_index << std::endl;
    auto next = doOneStep(prev, now, cur_index);

    // next is in current domain
    if(this->isPoint(next, cur_index)) {
        prev = now;
        now = next;
    }
    // next point is outside -> need to change domain.
    else {
        bool change_failed = true;
        for(chart_index new_index = 0; new_index < this->atlas_size; new_index++) {
            if (this->changePointIndex(prev, cur_index, new_index) and
                this->changePointIndex(now, cur_index, new_index) )
            {

                auto alt_next = doOneStep(prev, now, new_index);
                // Everything is good
                if(this->isPoint(alt_next, new_index)) {
                    prev = now;
                    now = alt_next;
                    cur_index = new_index;

                    change_failed = false;
                    break;
                }
            }
        }
        if(change_failed) {
            //std::cout << "Cur_index: " << cur_index << ". " << std::endl;
            //std::cout << now.to_str() << std::endl;
            //std::cout << next.to_str() << std::endl;
            // exit(1);
            // throw std::logic_error("AbstractRiemannianManifold<N>::doOneStepWithChange : failed to change point index. This situation must be impossible, please check configuration of manifold");
        }
    }
}


template
CUDA_F double AbstractRiemannianManifold<3>::integrateAlongPath(genPoint<3>,
                                                  Vec<3>, size_t,
                                                  const genFunction<3, double> &,
                                                  const weightFunction &,
                                                  double step) const;

template
CUDA_F Vec<3> AbstractRiemannianManifold<3>::integrateAlongPath(genPoint<3>,
                                                  Vec<3>, 
                                                  size_t,
                                                  const genFunction<3, Vec<3>> &,
                                                  const weightFunction &,
                                                  double step) const;


template class AbstractRiemannianManifold<3>;
template class AbstractRiemannianManifold<2>;
template class AbstractRiemannianManifold<1>;


