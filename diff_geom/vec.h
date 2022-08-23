#ifndef VEC_H
#define VEC_H

#include <array>
#include <iostream>
#include <math.h>

//template <std::size_t N>
//using vec = std::array<double, N>;

template <std::size_t N>
class vec : public std::array<double, N>
{
public:
    std::string to_str();

    vec<N> operator+ (const vec<N>& B);
    vec<N> operator* (const double& a); // We have right multiplication

    double norm2();
    double norm();

    vec<N> nomalize();
};

#endif // VEC_H
