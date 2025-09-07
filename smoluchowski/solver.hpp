#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "operator.hpp"

/*!
    \file solver.hpp
    \brief Main Smoluchowski Solver classes

    The file defines Smoluchowski Solver wich can solve
    forward linear and nonlinear problems, linear adjoint and
    inverse problems of source function reconstruction.
    File also defines abstract source function class which can be
    as a general function with time support in interval [t1;t2]
    as a function with divided variables with fix time dependence. 
*/

//! Abstract class of source function of particles
/*! Class describe source function v(x,t) where
 * if v(x,t)>0 thereis a concentration particle x arrive at the system at the
 * moment t;
 * if v(x,t)<0 thereis a concentration particle x left the system at the
 * moment t where x is from [0;H], t in [0;T].
 * We suppose that v(x,t)=0 if t not in [m1 * tau; m2 * tau].
 */
class SourceFunction {
protected:
    unsigned N, M;
    unsigned m1, m2;
    double h, tau;
public:
    //! Initialize source function v
    /*! Initialize source function v on [0;H]x[0;T] with uniform grid size NxM
     * which such that that v(x,t)=0 if t not in [m1 * tau; m2 * tau].
     */
    SourceFunction(unsigned _N, double _H, unsigned _M, double _T,
            unsigned _m1, unsigned _m2);
    //! The destructor
    virtual ~SourceFunction() { }
    //! Return concentration particle of size xi * h at moment ti * tau
    virtual double operator()(unsigned xi, unsigned ti) const = 0;
    //! Get square L2 norm of source function
    virtual double SquareNorm() const = 0;
    //! Update source function in gradient method
    virtual void Update(double *q, double alpha, double dzeta) = 0;
    //! Get size of uniform grid of particles
    unsigned GetN() const { return N; }
    //! Get size of uniform grid of time
    unsigned GetM() const { return M; }
    //! Get left bound of time support of source function
    unsigned Get_m1() const { return m1; }
    //! Get right bound of time support of source function
    unsigned Get_m2() const { return m2; }
    //! Get step of uniform grid of time
    double Get_tau() const { return tau; }
};

//! General source function of particles
/*! For details see SourceFunction. */
class SourceFunctionGeneral : public SourceFunction {
    double *v;
public:
    //! Initialize source function v
    /*! Initialize source function v on [0;H]x[0;T] with uniform grid size NxM
     * which such that that v(x,t)=0 if t not in [m1 * tau; m2 * tau] with
     * start value _v(x,t).
     * @param _N size of uniform grid on [0;H]
     * @param _H endpoint of particle segment
     * @param _M size of uniform grid on [0;T]
     * @param _T endpoint of time segment
     * @param _m1 left bound of time support of source function
     * @param _m2 right bound of time support of source function
     * @param _v double vector of size (N+1)*(M+1): _v[(N + 1) * k + i]=_v(i,k)
     */
    SourceFunctionGeneral(unsigned _N, double _H, unsigned _M, double _T,
            unsigned _m1, unsigned _m2, double *_v);
    //! The destructor
    ~SourceFunctionGeneral();
    double operator()(unsigned xi, unsigned ti) const;
    double SquareNorm() const;
    //! Update source function in gradient method
    /*! Update v <- v - dzeta * (alpha * v + q) */
    void Update(double *q, double alpha, double dzeta);
};

//! Source function with divided variables with fix time dependence
/*! Class describe source function v(x,t) with divided variables with fix
 * time dependence that is v(x,t) = vx(x) * vt(t) where vt is known.
 */
class SourceFunctionFixTime : public SourceFunction {
    double *vx;
    double *vt;
    double *vt_scale;
    double square_norm_vt;
public:
    //! Initialize source function v
    /*! Initialize values of source function
     * v(x,t) = vx(x) * vt(t) on [0;H]x[0;T]
     * with uniform grid size NxM which such that that vt(t)=0 if
     * t not in [m1 * tau; m2 * tau] with start value vx=_vx with fixed vt.
     * @param _N size of uniform grid on [0;H]
     * @param _H endpoint of particle segment
     * @param _M size of uniform grid on [0;T]
     * @param _T endpoint of time segment
     * @param _m1 left bound of time support of source function
     * @param _m2 right bound of time support of source function
     * @param _vt double vector of size M+1 of values _vt
     * @param _vx double vector of size N+1 of values _vx
     */
    SourceFunctionFixTime(unsigned _N, double _H, unsigned _M, double _T,
            unsigned _m1, unsigned _m2, double *_vt, double *_vx);
    //! The destructor
    ~SourceFunctionFixTime();
    double operator()(unsigned xi, unsigned ti) const;
    double SquareNorm() const;
    //! Update source function in gradient method
    /*! Update vx. Let t1 = m1 * tau, t2 = m2 * tau,
     * L2_inner(u,v) denote inner product in L2[t1;t2].
     * Then vx <- vx -
     * dzeta * (alpha * vx + L2_inner(q,vt) / L2_inner(vt,vt)).
     * */
    void Update(double *q, double alpha, double dzeta);
};


//! Class of stop criteria used in gradient method minimization of J_{alpha}(v)
//! TODO: rename this class
/*! Class describes stop criterias used in gradient method minimization of
 * Tikhonov regularization functional J_{alpha}(v) where
 * J_{alpha}(v) = alpha * square_norm_{L2([0;H]x[t1;t2])}(v) +
 * square_norm_{L2[0;H]}(qT) where v is source function,
 * qT(x)=h(x,T)-h_obs(x) where h is result
 * of solving forward linear problem.
 */
class StopCriteria {
public:
    //! StopCriteria mode
    enum stop_type {
        /*! stop by maximum number of iterations */
        stop_max_iters,
        /*! as stop_max_iters but also calculate minimized functional */
        stop_max_iters_calc_func,
        /*! stop by relative error of minimized functional */
        stop_func
    };
    //! Initialized stop criteria by max_iters
    /*! @param max_iter maximum number of iterations
     *  @param calc_func if true functional will be calculate on each iteration
     */
    StopCriteria(unsigned max_iter, bool calc_func = false);
    //! Initialized stop criteria by relative error of minimized functional
    /*! @param tol tolerance of relative error of functional
     */
    StopCriteria(double tol);
    //! Reset iterations to start
    void Reset();
    //! Get type of stop criteria
    stop_type GetType() const;
    //! Check is stop criteria satisfied
    bool IsStop() const;
    //! Update number of current iteration
    void UpdateCurIter() { cur_iter++; }
    //! Get number of current iteration
    unsigned GetCurIter() const { return cur_iter; }
    //! Check is calculating functional is needed
    bool NeedCalcFunc() const;
    //! Calculate functional J_{alpha}(v)
    /*! Calc J_{alpha}(v) = alpha * square_norm_{L2([0;H]x[t1;t2])}(v) +
     * square_norm_{L2[0;H]}(qT) where qT(x)=h(x,T)-h_obs(x) where h is result
     * of solving forward linear problem.
     * @param N size of uniform grid on particles
     * @param h step of uniform grid on particles
     * @param qT double vector of size N+1
     * @param alpha regularization parameter
     * @param square_norm_v equal square_norm_{L2([0;H]x[t1;t2])}(v) where v is
     * source function
     */
    void CalcFunc(unsigned N, double h, const double *qT,
        double alpha, double square_norm_v);
    //! Get value of functional on current iteration
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


//! Main Smoluchowski Solver class
/*! Class containes method to solve
 *  forward linear and nonlinear problems, linear adjoint and
 *  inverse problems of source function reconstruction.
 *  @note object of this class needeed only for solve inverse problem
 */
class SmoluchowskiSolver {
    SmoluchowskiLinearOperator& A;
    double *s, *h_obs;

    void calc_s();
public:
    //! Initialize solver of source function reconstruction problem
    //! by SmoluchowskiLinearOperator
    SmoluchowskiSolver(SmoluchowskiLinearOperator& _A);
    //! The destructor
    ~SmoluchowskiSolver();
    //! Get s(x) for details see description of SmoluchowskiLinearOperator
    const double *Get_s() const;
    //! Get size of uniform grid of particles
    unsigned GetN() const;

    //! Class of callback which can be used in inverse problem
    //! TODO: add callback to forward and adjoint problems
    class Callback {
    public:
        //! This method will be called on each iteration of inverse problem
        virtual void InverseProblemIter(const StopCriteria& stop_criteria) { }
    };
    
    //! Solve inverse problem of source function reconstruction
    /*! @param v initial value of source function
     *  @param g_obs vector of size N+1 of values difference between
     *  concentration of particles in start and end moments of time
     *  @param alpha regularization parameter (must be nonnegative)
     *  @param stopcriteria object of used stop criteria
     *  @param dzeta parameter of gradient method minimization of functional
     *  J_{alpha}(v) (for details see description of StopCriteria class)
     *  @param callback optional parameter to use information saved in
     *  stopcriteria object
     * @return v found reconstuction of source function
     */
    void inverse_problem(SourceFunction& v,
        const double *g_obs, double alpha, double dzeta,
        StopCriteria& stopcriteria, Callback* callback = 0);

    //! Abstract class of rightside which can be used in smoluchowki problems
    //! TODO: remove this class
    class RightSide {
    public:
        //! Get value of rightside in point (xi * h, ti * tau) where h,tau is
        //! steps of grids on particles and time
        virtual double operator()(unsigned xi, unsigned ti) const = 0;
    };

    //! Solve linear or nonlinear forward smoluchowski problem
    /*! Let A be a linear or nonlinear smoluchowski operator.
     * Method solve equation dg/dt + Ag = rightside with known g(t=0)
     * by Euler method.
     * @param A smoluchowski operator (linear or nonlinear)
     * @param T end time moment
     * @param M size of uniform grid on [0;T]
     * @param rightside right side of equation above
     * @param res vector of size (N+1)x(M+1) contains initial values of g:
     * g(i * h, 0) = res[i]
     * @return res vector of size (N+1)x(M+1) contains values of g on
     * grid [0;H]x[0;T]: res[(N + 1) * k + i] = g(i * h, k * tau).
     */
    static void forward_problem(SmoluchowskiOperator& A, double T, 
        unsigned M, const RightSide& rightside, double *res);

    //! Solve linear adjoint smoluchowski problem
    /*! Let A be a linear operator.
     * Method solve equation -dq/dt + A^{*}q = 0, with known g(t=T)
     * by Euler method.
     * @param A linear smoluchowski operator
     * @param T end time moment
     * @param M size of uniform grid on [0;T]
     * @param res vector of size (N+1)x(M+1) contains known g(t=T):
     * g(i * h, T) = res[(N+1) * M + i]
     * @return res vector of size (N+1)x(M+1) contains values of g on
     * grid [0;H]x[0;T]: res[(N + 1) * k + i] = g(i * h, k * tau).
     */
    static void adjoint_problem(SmoluchowskiLinearOperator& A, 
        double T, unsigned M, double *res);
};

//! Class of rightside without dependent of time which can be used in
//! smoluchowki problems
class RightSideNoTime : public SmoluchowskiSolver::RightSide {
    unsigned N;
    const double *arr; // no copy
public:
    //! Initialize rightside without dependence of time by vector of
    //! particles dependence
    /*! @param _N size of uniform particles grid
     * @param _arr double vector of size N+1 describes dependence righside by
     * particles
     * @note _arr saved without copy so don't change it while computations
     */
    RightSideNoTime(unsigned _N, const double *_arr);
    //! Set new righside
    void Set(unsigned _N, const double *_arr);
    double operator()(unsigned xi, unsigned ti) const;
};

//! Class of general rightside which can be used in smoluchowki problems
class RightSideGeneral : public SmoluchowskiSolver::RightSide {
    unsigned N, M;
    const double *arr; // no copy
public:
    //! Initialize general rightside
    /*! @param _N size of uniform particles grid
     * @param _M size of uniform time grid
     * @param _arr double vector of size (N+1)x(M+1) of values rightsize on
     * grid on [0;H]x[0;T]: rightside(xi * h,ti * tau)=arr[(N + 1) * ti + xi]
     * @note _arr saved without copy so don't change it while computations
     */
    RightSideGeneral(unsigned _N, unsigned _M, const double *_arr);
    //! Set new righside
    void Set(unsigned _N, unsigned _M, const double *_arr);
    double operator()(unsigned xi, unsigned ti) const;
};

//! Class of rightside generated by source function which can be used in
//! smoluchowki problems
class RightSideSourceFunction : public SmoluchowskiSolver::RightSide {
    const SourceFunction& v;
public:
    //! Initialize rightside by source function
    RightSideSourceFunction(const SourceFunction& _v);
    double operator()(unsigned xi, unsigned ti) const;
};

#endif /* SOLVER_HPP */
