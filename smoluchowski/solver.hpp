#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "operator.hpp"

class SourceFunction {
protected:
    unsigned N, M;
    unsigned m1, m2;
    double h, tau;
public:
    SourceFunction(unsigned _N, double _H, unsigned _M, double _T,
            unsigned _m1, unsigned _m2);
    virtual ~SourceFunction() { }
    virtual double operator()(unsigned xi, unsigned ti) const = 0;
    virtual double SquareNorm() const = 0;
    virtual void Update(double *q, double alpha, double dzeta) = 0;
    unsigned GetN() const { return N; }
    unsigned GetM() const { return M; }
    unsigned Get_m1() const { return m1; }
    unsigned Get_m2() const { return m2; }
    double Get_tau() const { return tau; }
};

class SourceFunctionGeneral : public SourceFunction {
    double *v;
public:
    SourceFunctionGeneral(unsigned _N, double _H, unsigned _M, double _T,
            unsigned _m1, unsigned _m2, double *_v);
    ~SourceFunctionGeneral();
    double operator()(unsigned xi, unsigned ti) const;
    double SquareNorm() const;
    void Update(double *q, double alpha, double dzeta);
};

class SourceFunctionFixTime : public SourceFunction {
    double *vx;
    double *vt;
    double *vt_scale;
    double square_norm_vt;
public:
    SourceFunctionFixTime(unsigned _N, double _H, unsigned _M, double _T,
            unsigned _m1, unsigned _m2, double *_vt, double *_vx);
    ~SourceFunctionFixTime();
    double operator()(unsigned xi, unsigned ti) const;
    double SquareNorm() const;
    void Update(double *q, double alpha, double dzeta);
};


class StopCriteria {
public:
    enum stop_type { stop_max_iters, stop_max_iters_calc_func, stop_func };
    StopCriteria(unsigned max_iter, bool calc_func = false);
    StopCriteria(double tol);
    void Reset();

    stop_type GetType() const;

    bool IsStop() const;

    void UpdateCurIter() { cur_iter++; }
    unsigned GetCurIter() const { return cur_iter; }

    bool NeedCalcFunc() const;
    void CalcFunc(unsigned N, double h, const double *qT,
        double alpha, double square_norm_v);
    double GetCurFunc() const { return func_new; }
private:
    char type;
    unsigned cur_iter;
    union un_param {
        unsigned max_iter;
        double tol;
    } param;

    double func_old, func_new;
};



/* create object only if you need solve inverse problem */
class SmoluchowskiSolver {
    SmoluchowskiLinearOperator& A;
    double *s, *h_obs;

    void calc_s();
public:
    SmoluchowskiSolver(SmoluchowskiLinearOperator& _A);
    ~SmoluchowskiSolver();
    const double *Get_s() const;
    unsigned GetN() const;

    class Callback {
    public:
        virtual void InverseProblemIter(const StopCriteria& stop_criteria) { }
    };

    void inverse_problem(SourceFunction& v,
        const double *g_obs, double alpha, double dzeta,
        StopCriteria& stopcriteria, Callback* callback = 0);

    //TODO: remove this class
    class RightSide {
    public:
        virtual double operator()(unsigned xi, unsigned ti) const = 0;
    };

    /* solve eq: dg/dt + Ag = rightside, g(t=0) = res(0 : N + 1) by Euler
    * g(t = k) = res((N + 1) * k : (N + 1) * (k + 1))
    * rightside, res are arrays with ((N + 1) * (M + 1)) elements  */
    static void forward_problem(SmoluchowskiOperator& A, double T, 
        unsigned M, const RightSide& rightside, double *res);

    /* solve eq by Euler: -dq/dt + A*q = 0, with
     * g(t=T) = res((N + 1) * M : (N + 1) * (M + 1));
     * q(t = k) = res((N + 1) * k : (N + 1) * (k + 1))
     * rightside, res are arrays with ((N + 1) * (M + 1)) elements  */
    static void adjoint_problem(SmoluchowskiLinearOperator& A, 
        double T, unsigned M, double *res);
};

class RightSideNoTime : public SmoluchowskiSolver::RightSide {
    unsigned N;
    double *arr; // no copy
public:
    RightSideNoTime(unsigned _N, double *_arr);
    void Set(unsigned _N, double *_arr);
    double operator()(unsigned xi, unsigned ti) const;
};

class RightSideGeneral : public SmoluchowskiSolver::RightSide {
    unsigned N, M;
    double *arr; // no copy
public:
    RightSideGeneral(unsigned _N, unsigned _M, double *_arr);
    void Set(unsigned _N, unsigned _M, double *_arr);
    double operator()(unsigned xi, unsigned ti) const;
};

class RightSideSourceFunction : public SmoluchowskiSolver::RightSide {
    const SourceFunction& v;
public:
    RightSideSourceFunction(const SourceFunction& _v);
    double operator()(unsigned xi, unsigned ti) const;
};

#endif /* SOLVER_HPP */
