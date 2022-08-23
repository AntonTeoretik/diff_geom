#include "vec.h"


template<std::size_t N>
std::string vec<N>::to_str()
{
    std::string res = "";
    for (auto c : (*this))
        res += std::to_string(c) + " ";
    return res;
}

template<std::size_t N>
vec<N> vec<N>::operator+(const vec<N> &B)
{
    vec<N> res;
    for(std::size_t i = 0; i < N; i++)
        res[i] = (*this)[i] + B[i];
    return res;
}

template<std::size_t N>
vec<N> vec<N>::operator*(const double &a)
{
    vec<N> res;
    for(std::size_t i = 0; i < N; i++)
        res[i] = a * (*this)[i];
    return res;
}

template<std::size_t N>
double vec<N>::norm2()
{
    double res = 0;
    for (auto c : (*this))
        res += c*c;
    return res;
}

template<std::size_t N>
double vec<N>::norm()
{
    return sqrt(norm2());
}

template<std::size_t N>
vec<N> vec<N>::nomalize()
{
    return (*this) * (1 / norm());
}


template class vec<3>;
template class vec<2>;
