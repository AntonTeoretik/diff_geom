#include "ellipsoid.h"
#include <iostream>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

using namespace std;

#include "algebra.h"
#include "metric.h"
#include "graphics.h"
#include "color.h"

#include <fstream>

#include <array>


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


class my_test_funtion_on_manifold : public genFunction<3, Vec<3>>
{
    Ellipsoid<3> S;

    CUDA_F double closeness_to_int(double x) const {
        double x1 = std::round(x);
        x1 = abs(x1 - x);
        return x1 <= 0.1 ? 2.0 : 0.0;
    }
    // genFunction interface
public:
    CUDA_F my_test_funtion_on_manifold(const Ellipsoid<3>& S) : S(S) {}

    CUDA_F Vec<3> operator ()(const genPoint<3> & p) const
    {
        Point<3> pp = p.p;

        if(p.i == 1 and !S.changePointIndex(pp, 1, 0)) {
            return Vec<3>::zero();
        }
        //pp.scale(0.1);

        // double x = closeness_to_int( pp[0] - 0.5);
        // double y = closeness_to_int( pp[1] - 0.5);
        // double z = closeness_to_int( pp[2] - 0.5);

        // return Vec<3>({x*y, y*z, z*x});

        pp.scale(4);

        int res = int(0.1 - abs(pp[0]) > 0) + int(0.1 - abs(pp[1]) > 0) + int(0.1 - abs(pp[2]) > 0);
        if (res < 1) {
            return {0.0, 0.0, 0.0};
        }
        double x = closeness_to_int(pp[0]);
        double y = closeness_to_int(pp[1]);
        double z = closeness_to_int(pp[2]);
        //double z = std::max(0.0, 0.1 - abs(pp[2]));

        return Vec<3>{x*y, y*z, z*x}*3.6;
    }
};

__global__ void render_kernel(double* out_img, int width, int height, int block_size) {
    int task_idx = threadIdx.x+blockIdx.x*blockDim.x;
    
    // Task indexing
    int nb_h = height / block_size;
    int nb_w = width / block_size;

    if (task_idx >= nb_h * nb_w) {
        return;
    }

    int yb = task_idx / nb_w;
    int xb = task_idx % nb_w;
    
    int x1 = xb * block_size;
    int x2 = (xb + 1) * block_size;

    int y1 = yb * block_size;
    int y2 = (yb + 1) * block_size;

    // Setup env
    Ellipsoid<3> S(Vec<4>({1.0, 1.0, 1.0, 1.0}), 9.0);
    my_test_funtion_on_manifold f(S);

    standartWeight w(0.2);

    //auto metr = EllipsoidMetric<3>({1,1,1,1}, NORTH);
    //std::cout << "Metr " << metr({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}) << std::endl ;

    Renderer rend(S, f, {0, {5,5,1}}, w);

    rend.number_of_points = 330;
    rend.step = 0.01;
    rend.screen.h_resolution = height;
    rend.screen.w_resolution = width;
    //1366x768

    rend.screen.height = 7.68;
    rend.screen.width = 13.66;
    rend.screen.distance_to_screen = 3.0;
    rend.setBasis(Vec<3>({-5, -5, -1}), Vec<3>({0, 0, -1}), Vec<3>({0, 1, 0}));

    // Do the work
    for (int x = x1; x < x2; ++x) {
        for (int y = y1; y < y2; ++y) {
            auto res = rend.trace_one_pixel(x, y);

            int offset = (y * width + x) * 3;
            out_img[offset] = res[0];
            out_img[offset + 1] = res[1];
            out_img[offset + 2] = res[2];
        }
    }
}

void test_one() {
    int width = 100;
    int height = 100;
    int y = 50;
    int x = 50;

    Ellipsoid<3> S(Vec<4>({1.0, 1.0, 1.0, 1.0}), 9.0);
    my_test_funtion_on_manifold f(S);
    standartWeight w;

    //auto metr = EllipsoidMetric<3>({1,1,1,1}, NORTH);
    //std::cout << "Metr " << metr({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}) << std::endl ;

    Renderer rend(S, f, {0, {0,0,0}}, w);

    rend.number_of_points = 500;
    rend.step = 0.01;
    rend.screen.h_resolution = width;
    rend.screen.w_resolution = height;
    rend.setBasis({1,3,2}, {3,2,1}, {1, 0, 0});

    auto res = rend.trace_one_pixel(x, y);
    printf("%f %f %f\n", res[0], res[1], res[2]);
}

int main()
{
    // test_one();
    // return 0;
    
    int width = 1280;
    int height = 720;
    
    thrust::device_vector<double> pixels(width * height * 3);
    gpuErrchk( cudaDeviceSynchronize() );

    // Set stack limit
    auto error = cudaDeviceSetLimit(cudaLimitStackSize, 16 * 1024);
    if (error != cudaSuccess) {
        printf("cudaDeviceSetLimit failed with %d, line(%d)\n", error, __LINE__);
        exit(EXIT_FAILURE);
    }

    size_t p_val;
    error = cudaDeviceGetLimit(&p_val, cudaLimitStackSize);
    if (error != cudaSuccess) {
        printf("cudaDeviceGetLimit failed with %d, line(%d)\n", error, __LINE__);
        exit(EXIT_FAILURE);
    }

    printf("stack size limit: %lu\n", p_val);
    
    int bs = 8;

    int block_dim = 32;
    int grid_dim = ((width / bs) * (height / bs)) / block_dim + 1;

    double* raw_gpu_pixels = thrust::raw_pointer_cast( pixels.data() );

    render_kernel<<<grid_dim, block_dim>>>(raw_gpu_pixels, width, height, bs);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

    thrust::host_vector<double> pixels_host = pixels;
    gpuErrchk( cudaDeviceSynchronize() );

    bitmap_image img(width, height);

    for(size_t x = 0; x < width; x++) {
        for(size_t y = 0; y < height; y++) {
            int offset = (y * width + x) * 3;
            Vec<3> v;
            v[0] = pixels_host[offset];
            v[1] = pixels_host[offset + 1];
            v[2] = pixels_host[offset + 2];
            img.set_pixel(x, y, Color(v).toRGB());
        }
    }
    img.save_image("img.bmp");


    return 0;
}
