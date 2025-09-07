#include "solver.hpp"
#include "math.h"

/*===========SourceFunction===============*/


SourceFunction::SourceFunction(unsigned _N, double _H, 
        unsigned _M, double _T, unsigned _m1, unsigned _m2)
    : N(_N), M(_M), m1(_m1), m2(_m2), h(_H / double(_N)), tau(_T / double(_M))
{ }

SourceFunctionGeneral::SourceFunctionGeneral(unsigned _N, double _H, 
    unsigned _M, double _T, unsigned _m1, unsigned _m2, double *_v)
    : SourceFunction(_N, _H, _M, _T, _m1, _m2)
{
    unsigned i, k;
    v = new double[(N + 1) * (M + 1)];
    for(k = 0; k <= M; ++k) {
        for(i = 0; i <= N; ++i) {
            if((k < m1) || (m2 < k)) {
                v[(N + 1) * k + i] = 0;
            } else {
                v[(N + 1) * k + i] = _v[(N + 1) * k + i];
            }
        }
    }
}

SourceFunctionGeneral::~SourceFunctionGeneral()
{
    delete [] v;
}

double SourceFunctionGeneral::operator()(unsigned xi, unsigned ti) const
{
#if HANDLEXCEPTION
    if((xi > N) || (ti > M)) {
        throw "assume xi <= N and ti <= M";
    }
#endif
    return v[(N + 1) * ti + xi];
}

double SourceFunctionGeneral::SquareNorm() const
{
    unsigned i, k;
    double res, tmp;
    res = 0;
    for(k = m1; k <= m2; ++k) {
        tmp = 0;
        for(i = 0; i <= N; ++i) {
            tmp += v[(N + 1) * k + i] * v[(N + 1) * k + i];
        }
        tmp -= 0.5 * (v[(N + 1) * k] + v[(N + 1) * k + N]);
        if((k == m1) || (k == m2)) {
            tmp *= 0.5;
        }
        res += tmp;
    }
    return tau * h * res;
}

void SourceFunctionGeneral::Update(double *q, double alpha, double dzeta)
{
    unsigned i, k;
    for(k = m1; k <= m2; ++k) {
        for(i = 0; i <= N; ++i) {
            v[(N + 1) * k + i] -= dzeta * (alpha * 
                v[(N + 1) * k + i] + q[(N + 1) * k + i]);
        }
    }
}

SourceFunctionFixTime::SourceFunctionFixTime(unsigned _N, double _H, 
        unsigned _M, double _T, unsigned _m1, unsigned _m2,
        double *_vt, double *_vx)
    : SourceFunction(_N, _H, _M, _T, _m1, _m2)
{
    unsigned i;
    vx = new double[N + 1];
    vt = new double[M + 1];
    vt_scale = new double[M + 1];
    for(i = 0; i <= N; ++i) {
        vx[i] = _vx[i];
    }
    square_norm_vt = 0;
    for(i = m1; i <= m2; ++i) {
        square_norm_vt += _vt[i] * _vt[i];
    }
    square_norm_vt = tau * (square_norm_vt -
            0.5 * (_vt[m1] * _vt[m1] + _vt[m2] * _vt[m2]));
    for(i = 0; i <= M; ++i) {
        if((i < m1) || (m2 < i)) {
            vt[i] = 0;
            vt_scale[i] = 0;
        } else {
            vt[i] = _vt[i];
            vt_scale[i] = _vt[i] / square_norm_vt;
        }
    }
}

SourceFunctionFixTime::~SourceFunctionFixTime()
{
    delete [] vx;
    delete [] vt;
    delete [] vt_scale;
}

double SourceFunctionFixTime::operator()(unsigned xi, unsigned ti) const
{
#if HANDLEXCEPTION
    if((xi > N) || (ti > M)) {
        throw "assume xi <= N and ti <= M";
    }
#endif
    return vx[xi] * vt[ti];
}

double SourceFunctionFixTime::SquareNorm() const
{
    unsigned i;
    double square_norm_vx;
    square_norm_vx = 0;
    for(i = 0; i <= N; ++i) {
        square_norm_vx += vx[i] * vx[i];
    }
    square_norm_vx = h * (square_norm_vx -
            0.5 * (vx[0] * vx[0] + vx[N] * vx[N]));
    return square_norm_vx * square_norm_vt;
}

void SourceFunctionFixTime::Update(double *q, double alpha, double dzeta)
{
    double tmp;
    unsigned i, k;
    for(i = 0; i <= N; ++i) {
        tmp = 0;
        for(k = m1 + 1; k < m2; ++k) {
            tmp += vt_scale[k] * q[(N + 1) * k + i];
        }
        tmp += 0.5 * (vt_scale[m1] * q[(N + 1) * m1 + i] +
                vt_scale[m2] * q[(N + 1) * m2 + i]);
        vx[i] -= dzeta * (alpha * vx[i] + tau * tmp);
    }
}

/*===========RightSide===============*/

RightSideNoTime::RightSideNoTime(unsigned _N, const double *_arr)
    : N(_N), arr(_arr)
{ }

void RightSideNoTime::Set(unsigned _N, const double *_arr)
{
    N = _N;
    arr = _arr;
}

double RightSideNoTime::operator()(unsigned xi, unsigned ti) const
{
#if HANDLEXCEPTION
    if(xi > N) {
        throw "assume xi <= N";
    }
#endif
    return arr[xi];
}

RightSideGeneral::RightSideGeneral(unsigned _N, unsigned _M, const double *_arr)
    : N(_N), M(_M), arr(_arr)
{ }

void RightSideGeneral::Set(unsigned _N, unsigned _M, const double *_arr)
{
    N = _N;
    M = _M;
    arr = _arr;
}

double RightSideGeneral::operator()(unsigned xi, unsigned ti) const
{
#if HANDLEXCEPTION
    if((xi > N) || (ti > M)) {
        throw "assume xi <= N and ti <= M";
    }
#endif
    return arr[(N + 1) * ti + xi];
}

RightSideSourceFunction::RightSideSourceFunction(const SourceFunction& _v)
    : v(_v)
{ }

double RightSideSourceFunction::operator()(unsigned xi, unsigned ti) const
{
    return v(xi, ti);
}

/*===========StopCriteria===============*/

// TODO: func_new not init by 1

StopCriteria::StopCriteria(unsigned max_iter, bool calc_func)
    : type(calc_func ? stop_max_iters_calc_func : stop_max_iters),
    cur_iter(0), func_old(0), func_new(1)
{
    param.max_iter = max_iter;
}

StopCriteria::StopCriteria(double tol)
    : type(stop_func), cur_iter(0), func_old(0), func_new(1)
{
    param.tol = tol;
}

void StopCriteria::Reset()
{
    cur_iter = 0;
    func_old = 0;
    func_new = 1;
}

StopCriteria::stop_type StopCriteria::GetType() const
{
    switch(type) {
    case stop_max_iters:
        return stop_max_iters;
    case stop_max_iters_calc_func:
        return stop_max_iters_calc_func;
    case stop_func:
        return stop_func;
    default:
        /* imposible */
        return stop_func;
    }
}

bool StopCriteria::NeedCalcFunc() const
{
    return type != stop_max_iters;
}

bool StopCriteria::IsStop() const
{
    if(type == stop_func) {
        return fabs(func_old - func_new) < param.tol * func_old;
    } else {
        return cur_iter >= param.max_iter;
    }
}

void StopCriteria::CalcFunc(unsigned N, double h, const double *qT,
        double alpha, double square_norm_v)
{
    unsigned i;
    double square_norm_qT = 0;
    for(i = 0; i <= N; ++i) {
        square_norm_qT += qT[i] * qT[i];
    }
    square_norm_qT = h * (square_norm_qT - 0.5 * (qT[0] + qT[N]));
    func_old = func_new;
    func_new = alpha * square_norm_v + square_norm_qT;
}


/*===========SmoluchowskiSolver===============*/

/* solve eq: dg/dt + Ag = rightside, g(t=0) = res(0 : N + 1) by Euler
 * g(t = k) = res((N + 1) * k : (N + 1) * (k + 1))
 * res is array with ((N + 1) * (M + 1)) elements  */
void SmoluchowskiSolver::forward_problem(SmoluchowskiOperator& A, 
        double T, unsigned M, const RightSide& rightside, double *res)
{
    unsigned i, k, N;
    double tau = T / double(M);
    N = A.GetSmolCalc().GetN();
    for(k = 0; k < M; ++k) {
        A.Apply(res + (N + 1) * k, res + (N + 1) * (k + 1));
        for(i = 0; i <= N; ++i) {
            res[(N + 1) * (k + 1) + i] = res[(N + 1) * k + i] + 
                tau * (rightside(i, k) - 
                        res[(N + 1) * (k + 1) + i]);
        }
    }
}

/* solve eq: -dq/dt + A*q = 0, g(t=T) = res((N + 1) * M : (N + 1) * (M + 1))
 * by Euler; q(t = k) = res((N + 1) * k : (N + 1) * (k + 1))
 * res is array with ((N + 1) * (M + 1)) elements  */
void SmoluchowskiSolver::adjoint_problem(SmoluchowskiLinearOperator& A, 
        double T, unsigned M, double *res)
{
    unsigned i, k, N;
    double tau = T / double(M);
    N = A.GetSmolCalc().GetN();
    for(k = M; k > 0; --k) {
        A.ApplyAdjoint(res + (N + 1) * k, res + (N + 1) * (k - 1));
        for(i = 0; i <= N; ++i) {
            res[(N + 1) * (k - 1) + i] = res[(N + 1) * k + i] -
                tau * res[(N + 1) * (k - 1) + i];
        }
    }
}


SmoluchowskiSolver::SmoluchowskiSolver(SmoluchowskiLinearOperator& _A)
    : A(_A)
{
    unsigned N = A.GetSmolCalc().GetN();
    h_obs = new double[N + 1];
    calc_s();
}

SmoluchowskiSolver::~SmoluchowskiSolver()
{
    delete [] s;
    delete [] h_obs;
}

void SmoluchowskiSolver::calc_s()
{
    unsigned i, N;
    const double *a, *c0;
    SmoluchowskiCalc& smol_calc = A.GetSmolCalc();
    a = A.Get_a();
    c0 = A.Get_c0();
    N = smol_calc.GetN();
    s = new double[N + 1];
    smol_calc.calc_L4(c0, s);
    smol_calc.calc_L1_fix_f(c0, h_obs);// A already fix c0; h_obs is tmp
    for(i = 0; i <= N; ++i) {
        s[i] += 0.5 * h_obs[i] - a[i] * c0[i];
    }
}

const double *SmoluchowskiSolver::Get_s() const
{
    return s;
}

unsigned SmoluchowskiSolver::GetN() const
{
    return A.GetSmolCalc().GetN();
}

void SmoluchowskiSolver::inverse_problem(SourceFunction& v,
        const double *g_obs, double alpha, double dzeta,
        StopCriteria& stop_criteria, Callback* callback)
{
    unsigned i, N, M;
    double T, h, *tmp_res;
    N = v.GetN();
    M = v.GetM();
    T = v.Get_tau() * M;
    h = A.Get_h();
    tmp_res = new double[(N + 1) * (M + 1)];
    for(i = 0; i <= N; ++i) {
        tmp_res[i] = 0;
    }
    /* calc g_1 */
    forward_problem(A, T, M, RightSideNoTime(N, s), tmp_res);
    /* calc h_obs = g_obs - g1(t = T) */
    for(i = 0; i <= N; ++i) {
        h_obs[i] = g_obs[i] - tmp_res[(N + 1) * M + i];
    }
    do {
        /* calc h */
        for(i = 0; i <= N; ++i) {
            tmp_res[i] = 0;
        }
        forward_problem(A, T, M, RightSideSourceFunction(v), tmp_res);
        /* calc q */
        for(i = 0; i <= N; ++i) {
            tmp_res[(N + 1) * M + i] -= h_obs[i];
        }
        if(stop_criteria.NeedCalcFunc()) {
            stop_criteria.CalcFunc(N, h, tmp_res + (N + 1) * M,
                    alpha, v.SquareNorm());
        }
        adjoint_problem(A, T, M, tmp_res);
        /* update v */
        v.Update(tmp_res, alpha, dzeta);
        if(callback) {
            callback->InverseProblemIter(stop_criteria);
        }
        stop_criteria.UpdateCurIter();
    } while (!stop_criteria.IsStop());
    delete [] tmp_res;
}
