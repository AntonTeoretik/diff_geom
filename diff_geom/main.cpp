#include <iostream>

using namespace std;

#include <algebra.h>
#include <geometry.h>


int main()
{
    Vec<3> A = {1,2,3};
    Vec<3> B = {1,2,3};
    cout << (A + B * 2).nomalized().to_str() << std::endl;

    MetricTensor<3> g;
    std::cout << g({0,0,0}, {0,1,0}, {0,1,1}) << std::endl;

    Matrix2D<3> M = { Vec<3>{0, 1, 1},
                      Vec<3>{1, 1, 8},
                      Vec<3>{0, 1, 3},
                    };

    std::cout << M.to_str() << std::endl;

    std::cout << M.inverse().to_str() << std::endl;

    std::cout << (M * M.inverse()).to_str() << std::endl;

    return 0;
}
