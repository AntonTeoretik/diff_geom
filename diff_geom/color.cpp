#include "color.h"

double Color::glow(double x)
{
    if (x < 0.0) return 0.0;
    return (1.0 - 1.0 / (x+1.0));
}

Color::Color()
{

}

Color::Color(Vec<3> v)
{
    v[0] = std::max(0.0, v[0]);
    v[1] = std::max(0.0, v[1]);
    v[2] = std::max(0.0, v[2]);

    double glow_red = glow(v[0] - 1.0);
    double glow_green = glow(v[1] - 1.0);
    double glow_blue = glow(v[2] - 1.0);

    r = std::min(v[0] + 0.5 * (glow_green + glow_blue), 1.0);
    g = std::min(v[1] + 0.5 * (glow_red + glow_blue), 1.0);
    b = std::min(v[2] + 0.5 * (glow_green + glow_red), 1.0);
}
