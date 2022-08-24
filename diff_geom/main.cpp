#include <iostream>

using namespace std;

#include <algebra.h>
#include <geometry.h>

#include <array>

int main()
{


    MetricTensor<2> g([](Point<2> P, Vec<2> v1, Vec<2> v2){
        Matrix2D<2> M = { {1, 0 },
                          {0, P[0]*P[0]} };

        return v1 * (M * v2);
    });

    cout << g.krist(1, 0, 1, {2,1}) << " " << endl;

    return 0;
}
