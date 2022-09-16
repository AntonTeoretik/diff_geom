#ifndef COLOR_H
#define COLOR_H

#include "algebra.h"

class Color
{
    double glow(double);

public:
    double r = 0.0, g = 0.0, b = 0.0;
    Color();
    Color(Vec<3>);

};

#endif // COLOR_H
