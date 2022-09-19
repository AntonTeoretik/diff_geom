#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <array>
#include <iostream>
#include <functional>
#include <math.h>

//template <std::size_t N>
//using vec = std::array<double, N>;

template <std::size_t N>
class Vec : public std::array<double, N>
{
public:
    std::string to_str() const;

    Vec<N> operator+ (const Vec<N>& B) const;
    Vec<N> operator- (const Vec<N>& B) const;

    Vec<N> operator* (const double& a) const; // We have right multiplication
    double operator* (const Vec<N>& B) const; // Scalar multiplication

    double norm2() const;
    double norm() const;

    Vec<N> normalized() const;

    void normalize();

    void addWithScale(const Vec<N>& v, double scale);
    void subtract(const Vec<N>& v);
    void scale(double a);

    static Vec<N> zero();
};

template <std::size_t N>
using Point = Vec<N>;

template <std::size_t N>
class Matrix2D : public std::array<std::array<double, N>, N>
{
public:
    //Matrix2D<N>();
    Matrix2D<N>() = default;
    Matrix2D<N>(std::initializer_list<std::array<double, N>>);
    Matrix2D<N>(std::function<double(int, int)> f);

    Matrix2D<N> inverse() const;
    Vec<N> operator* (const Vec<N>& A) const;
    Matrix2D<N> operator* (const Matrix2D<N>& M) const;
    Matrix2D<N> operator* (double a) const;


    std::string to_str() const;

};


template <std::size_t N>
std::array<Vec<N>, N> basis = []() {
    std::array<Vec<N>, N> res;
    for(std::size_t i = 0; i < N; i++) {
        for(std::size_t j = 0; j < N; j++) {
            res[i][j] = (i == j ? 1 : 0);
        }
    }
    return res;
}();




#endif // ALGEBRA_H
