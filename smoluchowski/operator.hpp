#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include "skeleton.hpp"
#include "fft.hpp"

class SmoluchowskiCalc {
protected:
    const Skeleton& K;
    const Skeleton& Psi;
    unsigned N;
    double H, h;
public:
    SmoluchowskiCalc(const Skeleton& K_skel, const Skeleton& Psi_skel, 
            double _H, unsigned _N);
    virtual ~SmoluchowskiCalc() { }
    unsigned GetN() const { return N; }
    double GetH() const { return H; }

    // int_{0}^{x}{K(x-y,y)f(x-y)g(y) dy}
    virtual void
    calc_L1(const double *f, const double *g, double *res) = 0;
    // int_{x}^{H}{K(y-x,x)f(y-x)g(y) dy}
    virtual void
    calc_L2(const double *f, const double *g, double *res) = 0;
    // save f (not copy)
    virtual void fix_f_arg(const double *f) = 0;
    virtual void calc_L1_fix_f(const double *g, double *res) = 0;
    virtual void calc_L2_fix_f(const double *g, double *res) = 0;

    // int_{0}^{x}{Psi(x,y)g(y) dy}
    virtual void calc_L3(const double *g, double *res) const = 0;
    // int_{x}^{H}{Psi(y,x)g(y) dy}
    virtual void calc_L4(const double *g, double *res) const = 0;
    // int_{0}^{H}{K(x,y)g(y) dy}
    virtual void calc_L5(const double *g, double *res) const = 0;
};

class SmoluchowskiCalcDirect : public SmoluchowskiCalc {
    double *Uf;
public:
    SmoluchowskiCalcDirect(const Skeleton& K_skel, 
            const Skeleton& Psi_skel, double _H, unsigned _N);
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

class SmoluchowskiCalcFast : public SmoluchowskiCalc {
    unsigned N_log2, N_2;
    FourierTransform Fourier;
    double *tmp_mv;   /* array with max(K.rank, Psi.rank) elem */
    double *Uf;
    Complex *Uf_fourier;
    Complex *tmp_complex0; // keep last part zero!
    Complex *tmp_complex1;
    Complex *tmp_complex2;
public:
    SmoluchowskiCalcFast(const Skeleton& K_skel, 
            const Skeleton& Psi_skel, double _H, unsigned _N);
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



class SmoluchowskiOperator {
protected:
    unsigned N;
    double H, h;
    SmoluchowskiCalc& smol_base;
    double *tmp_l;
public:
    SmoluchowskiOperator(SmoluchowskiCalc& _base);
    ~SmoluchowskiOperator();
    SmoluchowskiCalc& GetSmolCalc() { return smol_base; }
    double Get_h() const { return h; }
    virtual void Apply(const double *g, double *res) = 0;
};

class SmoluchowskiNonLinearOperator : public SmoluchowskiOperator {
public:
    SmoluchowskiNonLinearOperator(SmoluchowskiCalc& _base);
    ~SmoluchowskiNonLinearOperator();
    void Apply(const double *g, double *res);
};

class SmoluchowskiLinearOperator : public SmoluchowskiOperator {
    double *c0;
    double *a;
    
    void calc_a();
public:
    SmoluchowskiLinearOperator(SmoluchowskiCalc& _base, 
            const double *st); // with copy st
    ~SmoluchowskiLinearOperator();
    const double *Get_c0() const;
    const double *Get_a() const;
    virtual void Apply(const double *g, double *res);
    virtual void ApplyAdjoint(const double *g, double *res);
};

#endif /* OPERATOR_HPP */
