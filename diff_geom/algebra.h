#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <array>
#include <iostream>
#include <math.h>

//template <std::size_t N>
//using vec = std::array<double, N>;

template <std::size_t N>
class Vec : public std::array<double, N>
{
public:
    std::string to_str();

    Vec<N> operator+ (const Vec<N>& B);
    Vec<N> operator* (const double& a); // We have right multiplication

    double norm2();
    double norm();

    Vec<N> nomalize();
};

template <std::size_t N>
using Point = Vec<N>;

template <std::size_t N>
class Matrix2D : public std::array<std::array<double, N>, N>
{
public:
    Matrix2D<N> inverse();
    std::string to_str();

};

// Matrix inversion


#endif // ALGEBRA_H
