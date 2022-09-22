#ifndef WEIGHTS_H
#define WEIGHTS_H

#include <algorithm>

class weightFunction
{
public:
    virtual double operator () (double) const = 0;
};

class noWeight : public weightFunction
{
public:
    inline virtual double operator () (double) const {return 1.0;}
};

class linearWeight : public weightFunction
{
    double maxlen;
public:
    linearWeight(double maxlen) : maxlen(maxlen) {}
    inline virtual double operator () (double x) const {return std::max(1.0 - (x / maxlen), 0.0);}
};

class standartWeight : public weightFunction
{
    double scale;
public:
    standartWeight(double scale = 1.0) : scale(scale) {}
    inline virtual double operator () (double x) const {return 1.0 / (1.0 + scale * x);}
};


#endif // WEIGHTS_H
