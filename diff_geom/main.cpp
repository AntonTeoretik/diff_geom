#include <iostream>

using namespace std;

#include <algebra.h>
#include <geometry.h>


int main()
{

    auto f = [](Point<2> P){return P[0]*P[0] + P[1]*P[1];};

    cout << partialDer<2>({0.5,0}, basis<2>[0], f) << endl;

    return 0;
}
