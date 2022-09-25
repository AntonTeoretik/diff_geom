#ifndef WEIGHTS_H
#define WEIGHTS_H

#include <algorithm>
#include "util.h"

class weightFunction
{
public:
    virtual CUDA_F double operator () (double) const = 0;
};

class noWeight : public weightFunction
{
public:
    inline virtual CUDA_F double operator () (double) const {return 1.0;}
};

class linearWeight : public weightFunction
{
    double maxlen;
public:
    CUDA_F linearWeight(double maxlen) : maxlen(maxlen) {}
    inline virtual CUDA_F double operator () (double x) const {return fmax(1.0 - (x / maxlen), 0.0);}
};

class standartWeight : public weightFunction
{
    double scale;
public:
    CUDA_F standartWeight(double scale = 1.0) : scale(scale) {}
    inline virtual CUDA_F double operator () (double x) const {return 1.0 / (1.0 + scale * x);}
};


#endif // WEIGHTS_H
