#ifndef COLOR_H
#define COLOR_H

#include "algebra.h"
#include "bitimage/bitmap_image.hpp"

class Color
{
    static double glow(double x)
    {
        if (x < 0.0) return 0.0;
        return (1.0 - 1.0 / (x+1.0));
    }

public:
    double r = 0.0, g = 0.0, b = 0.0;
    Color();
    Color(Vec<3>);

    rgb_t toRB() const;

};

#endif // COLOR_H
