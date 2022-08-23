#include <iostream>

using namespace std;

#include <vec.h>


int main()
{
    vec<3> A = {1,2,3};
    vec<3> B = {1,2,3};
    cout << (A + B * 2).nomalize().to_str() << std::endl;



    return 0;
}
