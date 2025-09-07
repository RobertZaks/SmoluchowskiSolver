#ifndef FFT_HPP
#define FFT_HPP

/*!
    \file fft.hpp
    \brief Fourier Transform

    The file defines the Complex class and FourierTransform class
    with direct and fast matvec.
    See the FourierTransform class description for details.
*/

//! A usual class for a complex numbers
/*! Complex class implements a complex numbers as pair of two double values.
    Class contains a small number of operations required in program.
*/
class Complex {
    double re, im;
public:
    //! Make complex number from real and imaginary parts
    Complex(double _re = 0, double _im = 0);
    //! Get real part of complex number
    double GetRe() const;
    //! Get imaginary part of complex number
    double GetIm() const;
    //! Get squared (!) norm of complex number. TODO: rename method
    double GetNorm() const;
    //! Addition of complex numbers
    Complex operator+(const Complex &op2) const;
    //! Substraction of complex numbers
    Complex operator-(const Complex &op2) const;
    //! Multiplication of complex numbers
    Complex operator*(const Complex &op2) const;
    //! Division of complex numbers; check by yourself that op2.norm > 0
    Complex operator/(const Complex &op2) const;
    void operator+=(const Complex &op2);
    void operator*=(const Complex &op2);
    void operator/=(const Complex &op2);
};

//! A class of discrete Fourier Transform
/*! FourierTransform class implements a Cooley–Tukey algorithm for fast matvec
 *  with DFT matrix and direct matvec for testing.
*/
class FourierTransform {
    unsigned k;
    unsigned N;
    unsigned *rev;
    void set_reverse();
public:
    //! Initialize FourierTransform for DFT NxN matrix with N=2^k
    FourierTransform(unsigned _k);
    //! Type of matrix which will be use in operations: DFT or its conjugate
    enum mv_mode {normal, conj};
    //! Direct matvec DFT matrix with a vector
    /*!  @param vec Complex vector size N to witch DFT will apply
     *  @param res Complex vector size N where result will be save
     *  @param mode Calc matvec with DFT matrix or DFT conjugate
    */
    void Matvec(const Complex *vec, Complex *res, 
            mv_mode mode = normal) const;
    //! Fast matvec DFT matrix with a vector
    /*!  @param vec Complex vector size N to witch DFT will apply
     *  @param res Complex vector size N where result will be save
     *  @param mode Calc matvec with DFT matrix or DFT conjugate
    */
    void MatvecFast(const Complex *vec, Complex *res, 
            mv_mode mode = normal);
    //! The destructor
    ~FourierTransform();
};

#endif /* FFT_HPP */
