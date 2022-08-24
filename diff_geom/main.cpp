#include <iostream>

using namespace std;

#include <algebra.h>
#include <geometry.h>

#include <array>

int main()
{
    cout << Matrix2D<3>([](size_t i, size_t j){return i*i + j*j;}).to_str() << endl;

    return 0;
}
