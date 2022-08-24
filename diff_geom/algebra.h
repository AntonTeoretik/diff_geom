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
    std::string to_str() const;

    Vec<N> operator+ (const Vec<N>& B) const;
    Vec<N> operator* (const double& a) const; // We have right multiplication

    double norm2() const;
    double norm() const;

    Vec<N> nomalized() const;
};

template <std::size_t N>
using Point = Vec<N>;

template <std::size_t N>
class Matrix2D : public std::array<std::array<double, N>, N>
{
public:
    Matrix2D<N> inverse() const;
    Vec<N> operator* (const Vec<N>& A) const;
    Matrix2D<N> operator* (const Matrix2D<N>& M) const;
    Matrix2D<N> operator* (double a) const;


    std::string to_str() const;

};

// Matrix inversion


#endif // ALGEBRA_H
