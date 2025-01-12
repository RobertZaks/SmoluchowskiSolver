#include "fft.hpp"
#include <math.h>

/*=========Complex=======*/

Complex::Complex(double _re, double _im) 
    : re(_re), im(_im) { }

double Complex::GetRe() const
{
    return re;
}

double Complex::GetIm() const
{
    return im;
}

double Complex::GetNorm() const
{
    return re * re + im * im;
}

Complex Complex::operator+(const Complex &op2) const
{
    Complex res(re + op2.re, im + op2.im);
    return res;
}

Complex Complex::operator-(const Complex &op2) const
{
    Complex res(re - op2.re, im - op2.im);
    return res;
}

Complex Complex::operator*(const Complex &op2) const
{
    Complex res(re * op2.re - im * op2.im,
                re * op2.im + im * op2.re);
    return res;
}

Complex Complex::operator/(const Complex &op2) const
{
    double dvs = op2.GetNorm();
    Complex res((re * op2.re + im * op2.im) / dvs,
                (im * op2.re - re * op2.im) / dvs);
    return res;
}

void Complex::operator+=(const Complex &op2)
{
        re += op2.re;
        im += op2.im;
}

void Complex::operator*=(const Complex &op2)
{
    double tmp;
    tmp = re * op2.re - im * op2.im;
    im = re * op2.im + im * op2.re;
    re = tmp;
}

void Complex::operator/=(const Complex &op2) {
    double tmp, dvs = op2.GetNorm();
    tmp = (re * op2.re + im * op2.im) / dvs;
    im = (im * op2.re - re * op2.im) / dvs;
    re = tmp;
}

/*=====FourierTransform====*/

FourierTransform::FourierTransform(unsigned _k)
    : k(_k), rev(0)
{
    N = 1 << k;
}

FourierTransform::~FourierTransform()
{
    if(rev) {
        delete [] rev;
    }
}

void FourierTransform::set_reverse()
{
    unsigned i, h;
    h = -1;
    rev[0] = 0;
    for(i = 1; i < N; ++i) {
        /* h = floor(log_{2}{i}} */
        if((i & (i - 1)) == 0) {
            /* if i is a power of 2 */
            h++;
        }
        rev[i] = rev[i ^ (1 << h)];
        rev[i] |= (1 << (k - h - 1));
    }
}

void FourierTransform::Matvec(const Complex *vec, Complex *res, 
            mv_mode mode) const
{
    unsigned m, p;
    double a = -2.0 * M_PI / N;
    if(mode == conj) {
        a *= -1;
    }
    Complex w;
    for(m = 0; m < N; ++m) {
        res[m] = 0;
        for(p = 0; p < N; ++p) {
            w = Complex(cos(a * p * m), sin(a * p * m));
            res[m] += vec[p] * w;
        }
    }
}

void FourierTransform::MatvecFast(const Complex *vec, Complex *res, 
            mv_mode mode)
{
    unsigned m, p, j;
    Complex t, u, w, wm;
    double a = -2.0 * M_PI;
    if(mode == conj) {
        a *= -1;
    }
    if(!rev) {
        rev = new unsigned[N];
        set_reverse();
    }
    for(m = 0; m < N; ++m) {
        res[m] = vec[rev[m]];
    }
    
    for(m = 1; m <= N; m *= 2) {
        wm = Complex(cos(a / m), sin(a / m));
        for(p = 0; p < N; p += m) {
            w = 1;
            for(j = 0; j < m / 2; ++j) {
                t = w * res[p + j + m / 2];
                u = res[p + j];
                res[p + j] = u + t;
                res[p + j + m / 2] = u - t;
                w *= wm;
            }
        }
    }
}
