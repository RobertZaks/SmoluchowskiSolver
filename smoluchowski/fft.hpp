#ifndef FFT_HPP
#define FFT_HPP

class Complex {
    double re, im;
public:
    Complex(double _re = 0, double _im = 0);
    double GetRe() const;
    double GetIm() const;
    double GetNorm() const;
    Complex operator+(const Complex &op2) const;
    Complex operator-(const Complex &op2) const;
    Complex operator*(const Complex &op2) const;
    Complex operator/(const Complex &op2) const;
    void operator+=(const Complex &op2);
    void operator*=(const Complex &op2);
    void operator/=(const Complex &op2);
};

class FourierTransform {
    unsigned k;
    unsigned N;
    unsigned *rev;
    void set_reverse();
public:
    FourierTransform(unsigned _k);
    enum mv_mode {normal, conj};
    void Matvec(const Complex *vec, Complex *res, 
            mv_mode mode = normal) const;
    void MatvecFast(const Complex *vec, Complex *res, 
            mv_mode mode = normal);
    //void CalcDFT(const Complex *src, Complex *dst) const;
    //void CalcFFT(const Complex *src, Complex *dst);
    ~FourierTransform();
};

#endif /* FFT_HPP */
