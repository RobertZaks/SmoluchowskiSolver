#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include "skeleton.hpp"
#include "fft.hpp"

/*!
    \file operator.hpp
    \brief Smoluchowski operators

    The file defines abstract SmoluchowskiCalc class with methods required
    to apply smoluchowski operators and classes of smoluchowski operators.
*/


//! Abstract class of all calculation needed by smoluchowski operators applying
/*! Class defined functions to calculate results of applying
 *  operators Lk needed by smolychowski operators by trapezoidal rule.
 *  There are:
 *  - L1[f,g](x) = int_{0}^{x}{K(x-y,y)f(x-y)g(y)dy}
 *  - L2[f,g](x) = int_{x}^{H}{K(y-x,x)f(y-x)g(y)dy}
 *  - L3[g](x) = int_{0}^{x}{Psi(x,y)g(y)dy}
 *  - L4[g](x) = int_{0}^{H}{K(x,y)g(y)dy}
 *  - L5[g](x) = int_{0}^{H}{K(x,y)g(y)dy}
 */
class SmoluchowskiCalc {
protected:
    const Skeleton& K;
    const Skeleton& Psi;
    unsigned N;
    double H, h;
public:
    //! Initialize smoluchowski calculation
    /*! Let K, Psi be a coagulation and fragmentation kernel defined on
     * [0;H]x[0;H].
     * Let NxN be a size of uniform grid on [0;H]x[0;H].
     * @param K_skel skeleton approximation of matrix of K values on grid
     * @param Psi_skel skeleton approximation of matrix of Psi values on grid
     * @param _H endpoint of segment [0;H]
     * @param _N size of uniform grid on [0;H]
     */
    SmoluchowskiCalc(const Skeleton& K_skel, const Skeleton& Psi_skel, 
            double _H, unsigned _N);
    //! The destructor
    virtual ~SmoluchowskiCalc() { }
    //! Get the size of uniform grid
    unsigned GetN() const { return N; }
    //! Get the endpoint of segment [0;H]
    double GetH() const { return H; }

    //! Calculate int_{0}^{x}{K(x-y,y)f(x-y)g(y)dy}
    /*! Find res(x) = int_{0}^{x}{K(x-y,y)f(x-y)g(y)dy} by trapezoidal rule on
     *  uniform grid.
     *  @param f double vector values of function f of size N+1
     *  @param g double vector values of function g of size N+1
     *  @param res double vector values of function res of size N+1
     */
    virtual void
    calc_L1(const double *f, const double *g, double *res) = 0;
    //! Calculate int_{x}^{H}{K(y-x,x)f(y-x)g(y)dy}
    /*! Find res(x) = int_{x}^{H}{K(y-x,x)f(y-x)g(y)dy} by trapezoidal rule on
     *  uniform grid.
     *  @param f double vector values of function f of size N+1
     *  @param g double vector values of function g of size N+1
     *  @param res double vector values of function res of size N+1
     */
    virtual void
    calc_L2(const double *f, const double *g, double *res) = 0;
    //! Fix argument f in calc_L1 and calc_L2; see calc_L1_fix_f, calc_L2_fix_f
    virtual void fix_f_arg(const double *f) = 0;
    //! Calculate calc_L1 with previously saved argument f
    /*! Find res(x) = int_{0}^{x}{K(x-y,y)f(x-y)g(y)dy} by trapezoidal rule on
     *  uniform grid with f saved by fix_f_arg.
     *  @param g double vector values of function g of size N+1
     *  @param res double vector values of function res of size N+1
     */
    virtual void calc_L1_fix_f(const double *g, double *res) = 0;
    //! Calculate calc_L2 with previously saved argument f
    /*! Find res(x) = int_{x}^{H}{K(y-x,x)f(y-x)g(y)dy} by trapezoidal rule on
     *  uniform grid with f saved by fix_f_arg.
     *  @param g double vector values of function g of size N+1
     *  @param res double vector values of function res of size N+1
     */
    virtual void calc_L2_fix_f(const double *g, double *res) = 0;

    //! Calculate int_{0}^{x}{Psi(x,y)g(y)dy}
    /*! Find res(x) = int_{0}^{x}{Psi(x,y)g(y)dy} by trapezoidal rule on
     *  uniform grid.
     *  @param g double vector values of function g of size N+1
     *  @param res double vector values of function res of size N+1
     */
    virtual void calc_L3(const double *g, double *res) const = 0;
    //! Calculate int_{x}^{H}{Psi(y,x)g(y)dy}
    /*! Find res(x) = int_{x}^{H}{Psi(y,x)g(y)dy} by trapezoidal rule on
     *  uniform grid.
     *  @param g double vector values of function g of size N+1
     *  @param res double vector values of function res of size N+1
     */
    virtual void calc_L4(const double *g, double *res) const = 0;
    //! Calculate int_{0}^{H}{K(x,y)g(y)dy}
    /*! Find res(x) = int_{0}^{H}{K(x,y)g(y)dy} by trapezoidal rule on
     *  uniform grid.
     *  @param g double vector values of function g of size N+1
     *  @param res double vector values of function res of size N+1
     */
    virtual void calc_L5(const double *g, double *res) const = 0;
};

//! Class of direct calculation needed by smoluchowski operators applying
/*! Class contained computations of integral operators Lk
 *  (see SmoluchowskiCalc class for details)
 *  by directly trapezoidal rule.
 */
class SmoluchowskiCalcDirect : public SmoluchowskiCalc {
    double *Uf; // if fix_f_arg is called save K.GetU().elem(i,k) * f[i]
public:
    //! Initialize class of smoluchowski direct calculation
    /*! Let K, Psi be a coagulation and fragmentation kernel defined on
     * [0;H]x[0;H].
     * Let NxN be a size of uniform grid on [0;H]x[0;H].
     * @param K_skel skeleton approximation of matrix of K values on grid
     * @param Psi_skel skeleton approximation of matrix of Psi values on grid
     * @param _H endpoint of segment [0;H]
     * @param _N size of uniform grid on [0;H]
     */
    SmoluchowskiCalcDirect(const Skeleton& K_skel, 
            const Skeleton& Psi_skel, double _H, unsigned _N);
    //! The destructor
    ~SmoluchowskiCalcDirect();
    void calc_L1(const double *f, const double *g, double *res);
    void calc_L2(const double *f, const double *g, double *res);
    void fix_f_arg(const double *f);
    void calc_L1_fix_f(const double *g, double *res);
    void calc_L2_fix_f(const double *g, double *res);

    void calc_L3(const double *g, double *res) const;
    void calc_L4(const double *g, double *res) const;
    void calc_L5(const double *g, double *res) const;
};

//! Class of fast calculation needed by smoluchowski operators applying
/*! Class contained computations of integral operators Lk
 *  (see SmoluchowskiCalc class for details)
 *  by trapezoidal rule with fast algrorithm of matvec with matrix
 *  (or its transpose or their upper/lower triangular part) in skeleton format
 *  and fast discrete fourier transform.
 */
class SmoluchowskiCalcFast : public SmoluchowskiCalc {
    unsigned N_log2, N_2;
    FourierTransform Fourier;
    double *tmp_mv;   /* array with max(K.rank, Psi.rank) elem */
    // if fix_f_arg is called save K.GetU().elem(i,k) * f[i]
    double *Uf;
    // result of applying discrete fourier transform to Uf
    Complex *Uf_fourier;
    Complex *tmp_complex0; // keep last part zero!
    Complex *tmp_complex1;
    Complex *tmp_complex2;
public:
    //! Initialize class of smoluchowski fast calculation
    /*! Let K, Psi be a coagulation and fragmentation kernel defined on
     * [0;H]x[0;H].
     * Let NxN be a size of uniform grid on [0;H]x[0;H].
     * @param K_skel skeleton approximation of matrix of K values on grid
     * @param Psi_skel skeleton approximation of matrix of Psi values on grid
     * @param _H endpoint of segment [0;H]
     * @param _N size of uniform grid on [0;H]
     */
    SmoluchowskiCalcFast(const Skeleton& K_skel, 
            const Skeleton& Psi_skel, double _H, unsigned _N);
    //! The destructor
    ~SmoluchowskiCalcFast();
    void calc_L1(const double *f, const double *g, double *res);
    void calc_L2(const double *f, const double *g, double *res);
    void fix_f_arg(const double *f);
    void calc_L1_fix_f(const double *g, double *res);
    void calc_L2_fix_f(const double *g, double *res);

    void calc_L3(const double *g, double *res) const;
    void calc_L4(const double *g, double *res) const;
    void calc_L5(const double *g, double *res) const;
};


//! Abstract class of smoluchowski operators
class SmoluchowskiOperator {
protected:
    unsigned N;
    double H, h;
    SmoluchowskiCalc& smol_base;
    double *tmp_l;
public:
    //! Initialize smoluchowski operator by choose SmoluchowskiCalc realization
    SmoluchowskiOperator(SmoluchowskiCalc& _base);
    //! The destructor
    ~SmoluchowskiOperator();
    //! Get used SmoluchowskiCalc object
    SmoluchowskiCalc& GetSmolCalc() { return smol_base; }
    //! Get the step of uniform grid
    double Get_h() const { return h; }
    //! Apply operator
    /*! Applying smoluchowski operator to function g and result return in res
     * @param g double vector values of function g of size of uniform
     * grid (N+1)
     * @param res double vector values of function g of size of uniform
     * grid (N+1)
     */
    virtual void Apply(const double *g, double *res) = 0;
};

//! Nonlinear smoluchowski operator
/*! Let Lk be defined in SmoluchowskiCalc, then
 * SmoluchowskiNonLinearOperator.Apply(g)(x) = 
 * -1 * (0.5 * L1[g,g](x) - g(x) * L5[g](x) + L4[g](x) - g(x) / x * L3[y](x))
 * where L3[y] is L3 applying to function f(y)=y.
 */
class SmoluchowskiNonLinearOperator : public SmoluchowskiOperator {
    double *L3_y; // L3[y](x)
public:
    //! Initialize smoluchowski operator by choose SmoluchowskiCalc realization
    SmoluchowskiNonLinearOperator(SmoluchowskiCalc& _base);
    //! The destructor
    ~SmoluchowskiNonLinearOperator();
    void Apply(const double *g, double *res);
};

//! Linearized smoluchowski operator
/*! Let Lk be defined in SmoluchowskiCalc and L3[y]
 * is L3 applying to function f(y)=y.
 *  Let st be a distribution of particle in the moment of time t=0.
 *  Let a(x) = L5[st](x) + L3[y](x) / x,
 *  s(x) = 0.5 * L1[st,st](x) + L4[st](x) - a(x) * st(x).
 * Then SmoluchowskiLinearOperator.Apply(g)(x) =
 *  a(x) * g(x) - L1[st,g](x) + st(x) * L5[g](x) - L4[g](x)
 *  and SmoluchowskiLinearOperator.ApplyAdjoint(g)(x) =
 *  a(x) * g(x) + L5[st * g](x) - L2[st,g](x) - L3[g](x).
 */
class SmoluchowskiLinearOperator : public SmoluchowskiOperator {
    double *c0; // start distribution of particle
    double *a;
    
    void calc_a();
public:
    //! Initialize smoluchowski operator
    /*! @param _base choose SmoluchowskiCalc realization
     *  @param st set start distribution of particle
     */
    SmoluchowskiLinearOperator(SmoluchowskiCalc& _base, 
            const double *st); // with copy st
    //! The destructor
    ~SmoluchowskiLinearOperator();
    const double *Get_c0() const;
    const double *Get_a() const;
    virtual void Apply(const double *g, double *res);
    //! Apply adjoint to linearized smoluchowski operator
    /*! Applying adjoint to linearized smoluchowski operator to
     * function g and result return in res
     * @param g double vector values of function g of size of
     * uniform grid (N+1)
     * @param res double vector values of function g of size of
     * uniform grid (N+1)
     */
    virtual void ApplyAdjoint(const double *g, double *res);
};

#endif /* OPERATOR_HPP */
