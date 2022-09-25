#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <array>
#include <iostream>
// #include <functional>
#include <math.h>
#include "util.h"

//template <std::size_t N>
//using vec = std::array<double, N>;

template <std::size_t N>
class Vec 
{
public:
    std::string to_str() const;
    CUDA_F Vec<N>::Vec() {}
    CUDA_F Vec<N>::Vec(const double* d)  {
        for (int i = 0; i < N; ++i) {
            data[i] = d[i];
        }
    }

    CUDA_F Vec<N>::Vec(std::initializer_list<double> l)  {
        auto it = l.begin();
        for (int i = 0; i < N; ++i) {
            data[i] = *it;
            ++it;
        }
    }
    
    CUDA_F Vec<N> operator+ (const Vec<N>& B) const;
    CUDA_F Vec<N> operator- (const Vec<N>& B) const;

    CUDA_F Vec<N> operator* (const double& a) const; // We have right multiplication
    CUDA_F double operator* (const Vec<N>& B) const; // Scalar multiplication
    
    CUDA_F double& operator[] (size_t index) { return data[index]; }
    CUDA_F double operator[] (size_t index) const { return data[index]; }

    CUDA_F double norm2() const;
    CUDA_F double norm() const;

    CUDA_F Vec<N> normalized() const;

    CUDA_F void normalize();

    CUDA_F void addWithScale(const Vec<N>& v, double scale);
    CUDA_F void subtract(const Vec<N>& v);
    CUDA_F void scale(double a);

    CUDA_F static Vec<N> zero();
private:
    double data[N];
};

template <std::size_t N>
using Point = Vec<N>;

template <std::size_t N>
class Matrix2D 
{
public:
    //Matrix2D<N>();
    // Matrix2D<N>() = default;
    // Matrix2D<N>(std::initializer_list<std::array<double, N>>);
    // Matrix2D<N>(std::function<double(int, int)> f);
    using NNArray = double[N*N];

    CUDA_F Matrix2D<N>() {}
    // CUDA_F Matrix2D<N>(const NNArray* _data);
    
    template <typename Func>
    CUDA_F Matrix2D<N>(Func f) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i * N + j] = f(i, j);
            }
        }
    }

    CUDA_F Matrix2D<N> inverse() const;
    CUDA_F Vec<N> operator* (const Vec<N>& A) const;
    CUDA_F Matrix2D<N> operator* (const Matrix2D<N>& M) const;
    CUDA_F Matrix2D<N> operator* (double a) const;

    std::string to_str() const;
    CUDA_F double& get(int i, int j) { return data[i*N+j]; }
    CUDA_F const double& get(int i, int j) const { return data[i*N+j]; }

    double data[N*N];
};

template <std::size_t N>
inline CUDA_F void make_basis(Vec<N>* res) {
    for(std::size_t i = 0; i < N; i++) {
        for(std::size_t j = 0; j < N; j++) {
            res[i][j] = (i == j ? 1 : 0);
        }
    }
}

// template <std::size_t N>
// std::array<Vec<N>, N> basis = []() {
//     std::array<Vec<N>, N> res;
//     for(std::size_t i = 0; i < N; i++) {
//         for(std::size_t j = 0; j < N; j++) {
//             res[i][j] = (i == j ? 1 : 0);
//         }
//     }
//     return res;
// }();




#endif // ALGEBRA_H
