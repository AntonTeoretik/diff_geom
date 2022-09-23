#include "manifold.h"



template<std::size_t N>
std::vector<genPoint<N> > AbstractRiemannianManifold<N>::geodesic(
        const genPoint<N> &pt,
        Vec<N> dir,
        size_t num_of_pts,
        size_t dist,
        double step
        ) const
{
    Point<N> prev = pt.p;
    Point<N> now = prev + (dir * step); // Possible trouble here

    size_t count = 1;

    std::vector<genPoint<N> > res = {pt};
    if (dist == 1) {
        res.push_back({pt.i, now});
        count = 0;
    }

    chart_index cur_index = pt.i;

    for (size_t i = 0; i < num_of_pts; i++) {
        count++;
        doOneStepWithChange(prev, now, cur_index);

        if (count == dist) {
            count = 0;
            res.push_back({cur_index, now});
        }
    }
    return res;
}


template<std::size_t N>
std::vector<Vec<N> > AbstractRiemannianManifold<N>::orthogonalize(const genPoint<N> &pt, const std::vector<Vec<N> > &vecs, bool normalize) const
{
    return this->getMetric(pt.i).orthogonalize(pt.p, vecs, normalize);
}


template<class T, std::size_t N>
T integrateAlongPath(const std::vector<genPoint<N> > &points,
                     const std::function<T (genPoint<N>)> &func,
                     const std::function<double (double)> &weight,
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
        res = res + func(points[i]) * (step * weight(current_time));
    }
    return res;
}

template<size_t N>
template<class T>
T AbstractRiemannianManifold<N>::integrateAlongPath(genPoint<N> start,
                                                  Vec<N> dir, size_t num_of_pts,
                                                  const genFunction<N, T> &func,
                                                  const weightFunction &weight,
                                                  double step) const
{
    //std::cout << "AbstractRiemannianManifold<N>::integrateAlongPath : " << dir.to_str() << std::endl;

    double current_time = 0.0;

    chart_index cur_index = start.i;
    Point<N> prev = start.p;

    //NORMALIZE DIR
    double norm2 = getMetric(cur_index)(prev, dir, dir);
    dir.scale(1.0 / sqrt(norm2));

    Point<N> now = prev + (dir * step); // Possible trouble here

    T res = func(genPoint<N>{cur_index, now}) * weight(0.0) * step;

    for(size_t i = 1; i <= num_of_pts; i++) {
        current_time += step;
        doOneStepWithChange(prev, now, cur_index);
        res = res + func({cur_index, now}) * (step * weight(current_time));
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
void AbstractRiemannianManifold<N>::doOneStepWithChange(Point<N> &prev, Point<N> &now, chart_index &cur_index) const
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
            throw std::logic_error("AbstractRiemannianManifold<N>::doOneStepWithChange : failed to change point index. This situation must be impossible, please check configuration of manifold");
        }
    }
}


template
double AbstractRiemannianManifold<3>::integrateAlongPath(genPoint<3>,
                                                  Vec<3>, size_t,
                                                  const genFunction<3, double> &,
                                                  const weightFunction &,
                                                  double step) const;

template
Vec<3> AbstractRiemannianManifold<3>::integrateAlongPath(genPoint<3>,
                                                  Vec<3>, size_t,
                                                  const genFunction<3, Vec<3>> &,
                                                  const weightFunction &,
                                                  double step) const;


template class AbstractRiemannianManifold<3>;
template class AbstractRiemannianManifold<2>;
template class AbstractRiemannianManifold<1>;


