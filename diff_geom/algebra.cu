#include "algebra.h"


template<std::size_t N>
std::string Vec<N>::to_str() const
{
    std::string res = "";
    for (int i = 0; i < N; ++i)
        res += std::to_string(data[i]) + " ";
    return res;
}

template<std::size_t N>
CUDA_F Vec<N> Vec<N>::operator+(const Vec<N> &B) const
{
    Vec<N> res;
    for(std::size_t i = 0; i < N; i++)
        res[i] = data[i] + B[i];
    return res;
}

template<std::size_t N>
CUDA_F Vec<N> Vec<N>::operator-(const Vec<N> &B) const
{
    Vec<N> res;
    for(std::size_t i = 0; i < N; i++)
        res[i] = data[i] - B[i];
    return res;
}

template<std::size_t N>
CUDA_F Vec<N> Vec<N>::operator*(const double &a) const
{
    Vec<N> res;
    for(std::size_t i = 0; i < N; i++)
        res[i] = a * data[i];
    return res;
}

template<std::size_t N>
CUDA_F double Vec<N>::operator*(const Vec<N> &v) const
{
    double res = 0;
    for(std::size_t i = 0; i < N; i++)
        res += data[i] * v[i];
    return res;
}

template<std::size_t N>
CUDA_F double Vec<N>::norm2() const
{
    double res = 0;
    for (int i = 0; i < N; ++i) {
        double c = data[i];
        res += c*c;
    }
    return res;
}

template<std::size_t N>
CUDA_F double Vec<N>::norm() const
{
    return sqrt(norm2());
}

template<std::size_t N>
CUDA_F Vec<N> Vec<N>::normalized() const
{
    return (*this) * (1 / norm());
}

template<std::size_t N>
CUDA_F void Vec<N>::normalize()
{
    double inv_norm = 1.0 / this->norm();
    for(std::size_t i = 0; i < N; i++)
        data[i] *= inv_norm;
}

template<std::size_t N>
CUDA_F void Vec<N>::addWithScale(const Vec<N> &v, double scale)
{
    for(std::size_t i = 0; i < N; i++)
        data[i] += v[i] * scale;
}

template<std::size_t N>
CUDA_F void Vec<N>::subtract(const Vec<N> &v)
{
    for(std::size_t i = 0; i < N; i++)
        data[i] -= v[i];
}

template<std::size_t N>
CUDA_F void Vec<N>::scale(double a)
{
    for(std::size_t i = 0; i < N; i++)
        data[i] *= a;
}

template<std::size_t N>
CUDA_F Vec<N> Vec<N>::zero()
{
    Vec<N> res;
    for(std::size_t i = 0; i < N; i++)
        res[i] = 0;
    return res;
}



template<>
CUDA_F Matrix2D<1> Matrix2D<1>::inverse() const
{
    Matrix2D<1> mat;
    mat.get(0,0) = 1.0 / get(0,0);
    return mat;
}

template<>
CUDA_F Matrix2D<2> Matrix2D<2>::inverse() const
{
    double a = data[0];
    double b = data[1];
    double c = data[2];
    double d = data[3];

    double detinv = 1.0 / (a * d - b * c);

    Matrix2D<2> res;
    res.get(0,0) = d * detinv;
    res.get(0,1) = -b * detinv;
    res.get(1,0) = -c * detinv;
    res.get(1,1) = a * detinv;
    return res;
}


// template<std::size_t N>
// CUDA_F Matrix2D<N>::Matrix2D(std::initializer_list<std::array<double, N> > list)
// {
//     for(std::size_t i = 0; i < N; i++) {
//         (*this)[i] = std::array<double, N>{*(list.begin() + i)};
//     }
// }

// template<std::size_t N>
// Matrix2D<N>::Matrix2D(std::function<double (int, int)> f)
// {
//     for(std::size_t i = 0; i < N; i++) {
//         for(std::size_t j = 0; j < N; j++) {
//             (*this)[i][j] = f(i, j);
//         }
//     }
// }

// template <size_t N>
// CUDA_F Matrix2D<N>::Matrix2D(const NNArray* _data) {
//     for (int i = 0; i < N; ++i) {
//         for (int j = 0; j < N; ++j) {
//             get(i, j) = (*_data)[i][j];
//         }
//     }
// }

template<>
CUDA_F Matrix2D<3> Matrix2D<3>::inverse() const
{
    double a11 = data[0];
    double a12 = data[1];
    double a13 = data[2];

    double a21 = data[3];
    double a22 = data[4];
    double a23 = data[5];

    double a31 = data[6];
    double a32 = data[7];
    double a33 = data[8];

    double detinv = 1.0 / (a11*a22*a33 + a13*a21*a32 + a12*a23*a31 - a13*a22*a31 - a12*a21*a33 - a11*a23*a32);

    Matrix2D<3> adj;
    double inv_data[9] = {
        a22*a33 - a23*a32, -a12*a33 + a13*a32, a23*a12 - a22*a13,
        -a21*a33 + a31*a23,  a11*a33 - a31*a13, -a11*a23 + a21*a13,
        a21*a32 - a22*a31, -a11*a32 + a12*a31,  a22*a11 - a21*a12 
    };
    
    memcpy(adj.data, inv_data, 9 * sizeof(double));
    return adj * detinv;
}

template<std::size_t N>
CUDA_F Vec<N> Matrix2D<N>::operator*(const Vec<N> &A) const
{
    Vec<N> res;
    for(std::size_t i = 0; i < N; i++) {
        double subres = 0.0;
        for(std::size_t j = 0; j < N; j++)
            subres += get(i, j) * A[j];
        res[i] = subres;
    }
    return res;
}

template<std::size_t N>
CUDA_F Matrix2D<N> Matrix2D<N>::operator*(const Matrix2D<N> &A) const
{
    Matrix2D<N> res;
    for(std::size_t i = 0; i < N; i++) {
        for(std::size_t j = 0; j < N; j++) {
            double subres = 0.0;
            for(std::size_t r = 0; r < N; r++) {
                subres += get(i, r) * A.get(r, j);
            }
            res.get(i,j) = subres;
        }
    }

    return res;
}

template<std::size_t N>
CUDA_F Matrix2D<N> Matrix2D<N>::operator*(double a) const
{
    Matrix2D<N> res;
    for(std::size_t i = 0; i < N; i++) {
        for(std::size_t j = 0; j < N; j++) {
            res.get(i, j) = a * get(i, j);
        }
    }

    return res;
}

template<std::size_t N>
std::string Matrix2D<N>::to_str() const
{
    std::string str = "";
    for(std::size_t i = 0; i < N; i++) {
        for(std::size_t j = 0; j < N; j++) {
            str += std::to_string(get(i, j)) + " ";
        }
        str += "\n";
    }
    return str;
}



template class Vec<1>;
template class Vec<2>;
template class Vec<3>;
template class Vec<4>;

template class Matrix2D<1>;
template class Matrix2D<2>;
template class Matrix2D<3>;
template class Matrix2D<4>;



