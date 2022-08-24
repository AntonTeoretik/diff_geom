#include "algebra.h"


template<std::size_t N>
std::string Vec<N>::to_str() const
{
    std::string res = "";
    for (auto c : (*this))
        res += std::to_string(c) + " ";
    return res;
}

template<std::size_t N>
Vec<N> Vec<N>::operator+(const Vec<N> &B) const
{
    Vec<N> res;
    for(std::size_t i = 0; i < N; i++)
        res[i] = (*this)[i] + B[i];
    return res;
}

template<std::size_t N>
Vec<N> Vec<N>::operator*(const double &a) const
{
    Vec<N> res;
    for(std::size_t i = 0; i < N; i++)
        res[i] = a * (*this)[i];
    return res;
}

template<std::size_t N>
double Vec<N>::norm2() const
{
    double res = 0;
    for (auto c : (*this))
        res += c*c;
    return res;
}

template<std::size_t N>
double Vec<N>::norm() const
{
    return sqrt(norm2());
}

template<std::size_t N>
Vec<N> Vec<N>::nomalized() const
{
    return (*this) * (1 / norm());
}



template<>
Matrix2D<1> Matrix2D<1>::inverse() const
{
    return {Vec<1>{1.0 / (*this)[0][0]}};
}

template<>
Matrix2D<2> Matrix2D<2>::inverse() const
{
    double a = (*this)[0][0];
    double b = (*this)[0][1];
    double c = (*this)[1][0];
    double d = (*this)[1][1];

    double detinv = 1.0 / (a * d - b * c);

    return {
        Vec<2>{d * detinv, -b * detinv},
        Vec<2>{-c * detinv, a * detinv}
    };
}

//template<std::size_t N>
//Matrix2D<N>::Matrix2D<N>(std::function<double (int, int)> f)
//{

//}


template<std::size_t N>
Matrix2D<N>::Matrix2D(std::initializer_list<std::array<double, N> > list)
{
    for(std::size_t i = 0; i < N; i++) {
        (*this)[i] = std::array<double, N>{*(list.begin() + i)};
    }
}

template<std::size_t N>
Matrix2D<N>::Matrix2D(std::function<double (int, int)> f)
{
    for(std::size_t i = 0; i < N; i++) {
        for(std::size_t j = 0; j < N; j++) {
            (*this)[i][j] = f(i, j);
        }
    }
}

template<>
Matrix2D<3> Matrix2D<3>::inverse() const
{
    double a11 = (*this)[0][0];
    double a12 = (*this)[0][1];
    double a13 = (*this)[0][2];

    double a21 = (*this)[1][0];
    double a22 = (*this)[1][1];
    double a23 = (*this)[1][2];

    double a31 = (*this)[2][0];
    double a32 = (*this)[2][1];
    double a33 = (*this)[2][2];

    double detinv = 1.0 / (a11*a22*a33 + a13*a21*a32 + a12*a23*a31 - a13*a22*a31 - a12*a21*a33 - a11*a23*a32);

    Matrix2D<3> adj = {
        { a22*a33 - a23*a32, -a12*a33 + a13*a32, a23*a12 - a22*a13 },
        {-a21*a33 + a31*a23,  a11*a33 - a31*a13, -a11*a23 + a21*a13 },
        { a21*a32 - a22*a31, -a11*a32 + a12*a31,  a22*a11 - a21*a12 }
    };

    return adj * detinv;
}

template<std::size_t N>
Vec<N> Matrix2D<N>::operator*(const Vec<N> &A) const
{
    Vec<N> res;
    for(std::size_t i = 0; i < N; i++) {
        double subres = 0.0;
        for(std::size_t j = 0; j < N; j++)
            subres += (*this)[i][j] * A[j];
        res[i] = subres;
    }
    return res;
}

template<std::size_t N>
Matrix2D<N> Matrix2D<N>::operator*(const Matrix2D<N> &A) const
{
    Matrix2D<N> res;
    for(std::size_t i = 0; i < N; i++) {
        for(std::size_t j = 0; j < N; j++) {
            double subres = 0.0;
            for(std::size_t r = 0; r < N; r++) {
                subres += (*this)[i][r] * A[r][j];
            }
            res[i][j] = subres;
        }
    }
    return res;
}

template<std::size_t N>
Matrix2D<N> Matrix2D<N>::operator*(double a) const
{
    Matrix2D<N> res;
    for(std::size_t i = 0; i < N; i++) {
        for(std::size_t j = 0; j < N; j++) {
            res[i][j] = (*this)[i][j] * a;
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
            str += std::to_string((*this)[i][j]) + " ";
        }
        str += "\n";
    }
    return str;
}



template class Vec<2>;
template class Vec<3>;

template class Matrix2D<2>;
template class Matrix2D<3>;

