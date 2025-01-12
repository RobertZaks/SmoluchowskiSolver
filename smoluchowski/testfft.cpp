#include "fft.hpp"
#include <stdio.h>
#include <time.h>
#include <math.h>

enum { max_N = 4096 };

int main() {
    unsigned i, k, N;
    double x;
    Complex *vec, *res_fft, *res_dft, *res_ifft, *res_idft, tmp;
    clock_t stime;
    vec = new Complex[max_N];
    res_fft = new Complex[max_N];
    res_dft = new Complex[max_N];
    res_ifft = new Complex[max_N];
    res_idft = new Complex[max_N];
    
    srand(time(NULL));
    for(i = 0; i < max_N; ++i) {
        //scanf("%lf", &x);
        x = 10000.0 * rand() / (RAND_MAX + 1.0);
        vec[i] = x;
    }

    for(k = 0; (unsigned)(1 << k) <= max_N; ++k) {
        N = 1 << k;
        printf("---> N = %d\n", N);
        FourierTransform FT(k);
        stime = clock();
        FT.Matvec(vec, res_dft);
        printf("dft time:%ld\n", clock() - stime);
        stime = clock();
        FT.MatvecFast(vec, res_fft);
        printf("fft time:%ld\n", clock() - stime);
        x = 0;
        for(i = 0; i < N; ++i) {
            //printf("dft: %f %f\n", res_dft[i].GetRe(), res_dft[i].GetIm());
            //printf("fft: %f %f\n", res_fft[i].GetRe(), res_fft[i].GetIm());
            tmp = res_fft[i] - res_dft[i];
            x += tmp.GetNorm();
        }
        printf("||dft - fft||^2: %f\n", x);
        
        FT.Matvec(res_dft, res_idft, FourierTransform::conj);
        for(i = 0; i < N; ++i) {
            tmp = res_idft[i] / N- vec[i];
            x += tmp.GetNorm();
        }
        printf("||idft(dft(vec)) - vec||^2: %f\n", x);
        
        FT.MatvecFast(res_fft, res_ifft, FourierTransform::conj);
        for(i = 0; i < N; ++i) {
            tmp = res_ifft[i] / N - vec[i];
            x += tmp.GetNorm();
        }
        printf("||ifft(fft(vec)) - vec||^2: %f\n", x);
    }
    delete [] res_dft;
    delete [] res_fft;
    delete [] vec;
}
