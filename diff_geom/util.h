#pragma once

#define CUDA_F __host__ __device__


template <typename T>
struct Span {
    T* data;
    size_t size;
    
    CUDA_F T& operator[] (size_t i) {
        return data[i];
    }

    CUDA_F const T& operator[] (size_t i) const {
        return data[i];
    }
};