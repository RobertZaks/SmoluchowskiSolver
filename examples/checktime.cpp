#include <smoluchowski/matrix.hpp>
#include <smoluchowski/cross.hpp>
#include <smoluchowski/solver.hpp>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>

double K(unsigned i, unsigned j, double h)
{
#if 1
    return pow(pow((i + 1) * h, 1.0 / 3.0) + pow((j + 1) * h, 1.0 / 3.0), 2) *
        sqrt(1.0 / ((i + 1) * h) + 1.0 / ((j + 1) * h));
#endif
    //return 1.0;
    return pow((i + 1) * h * (j + 1) * h, 1.0 / 6.0);
}

double Psi(unsigned i, unsigned j, double h)
{
    //return 0.0;
    return exp(-3.0 * h * double(i));
}

double c0(double x)
{
    return exp(-x);
}

double v_solution(double x, double t)
{
    return exp(-x);
    return exp(-2.0 * x);
}

double v0(double x, double t)
{
    return exp(-2.0 * x);
    return exp(-x);
}

double time_diff(timeval &end, timeval &start)
{
    return (end.tv_sec - start.tv_sec) + 
        (end.tv_usec - start.tv_usec) * 1e-6;
}

class SmoluchowskiLinearOperatorSaveTime : public SmoluchowskiLinearOperator {
    struct timeval tv1, tv2;
    double max_apply;
    double max_apply_adjoint;
public:
    SmoluchowskiLinearOperatorSaveTime(SmoluchowskiCalc& _base,
            const double *st)
        : SmoluchowskiLinearOperator(_base, st), max_apply(0),
        max_apply_adjoint(0)
    { }
    ~SmoluchowskiLinearOperatorSaveTime()
    { }
    double GetMaxApplyTime() const { return max_apply; }
    double GetMaxApplyAdjointTime() const { return max_apply_adjoint; }
    void Apply(const double *g, double *res) {
        double tmp;
        gettimeofday(&tv1, 0);
        SmoluchowskiLinearOperator::Apply(g, res);
        gettimeofday(&tv2, 0);
        tmp = time_diff(tv2, tv1);
        if(tmp > max_apply) {
            max_apply = tmp;
        }
    }
    void ApplyAdjoint(const double *g, double *res) {
        double tmp;
        gettimeofday(&tv1, 0);
        SmoluchowskiLinearOperator::ApplyAdjoint(g, res);
        gettimeofday(&tv2, 0);
        tmp = time_diff(tv2, tv1);
        if(tmp > max_apply_adjoint) {
            max_apply_adjoint = tmp;
        }
    }
};

class CallbackPrintIter : public SmoluchowskiSolver::Callback {
public:
    CallbackPrintIter() { }
    ~CallbackPrintIter() { }
    void InverseProblemIter(const StopCriteria& stop_criteria) {
        fprintf(stderr, "iter: %d\n", stop_criteria.GetCurIter());
    }
};

int main()
{
    unsigned i, k;
    unsigned N = 1000, M = 80;
    unsigned m1, m2;
    unsigned num_iter;
    double H = 6.0, T = 0.5;
    double t1, t2, tau, h, alpha, c1, dzeta, tol_cross = 1e-6;
    double *res, *rightside, *ptmp;
    FILE *f;
    h = H / double(N);
    tau = T / double(M);
    rightside = new double[(N + 1) * (M + 1)];
    res = new double[(N + 1) * (M + 1)];
    ptmp = new double[N + 1];
    
    num_iter = 2;

    alpha = 0.05;
    c1 = T;
    dzeta = 2.0 / (2.0 * alpha + c1);
    t1 = T / 5;
    t2 = 3 * T / 5;
    m1 = floor(t1 / T * M);
    m2 = ceil(t2 / T * M);
    m1 = floor(t1 / T * M);
    m2 = ceil(t2 / T * M);
    printf("N: %d, H: %f, M: %d, T: %f\n", N, H, M, T);
    printf("m1=%d, m2=%d\n", m1, m2);

    FunctionMatrix K_M(N + 1, N + 1, h, K);
    Cross K_approx(K_M);
    K_approx.Approximate(tol_cross);
    FunctionMatrix Psi_M(N + 1, N + 1, h, Psi);
    Cross Psi_approx(Psi_M);
    Psi_approx.Approximate(tol_cross);
    printf("K   rank: %d\n", K_approx.GetRank());
    printf("Psi rank: %d\n", Psi_approx.GetRank());

    for(i = 0; i <= N; ++i) {
        ptmp[i] = c0(i * h);
        res[i] = 0;
    }
    
#ifdef CALC_DIRECT
    SmoluchowskiCalcDirect smol_base(K_approx, Psi_approx, H, N);
#else
    SmoluchowskiCalcFast smol_base(K_approx, Psi_approx, H, N);
#endif
    SmoluchowskiLinearOperatorSaveTime smol_op(smol_base, ptmp);
    SmoluchowskiSolver solver(smol_op);
    

    for(i = 0; i <= N; ++i) {
        for(k = 0; k <= M; ++k) {
           rightside[(N + 1) * k + i] = solver.Get_s()[i];
           if((m1 <= k) && (k <= m2)) {
                rightside[(N + 1) * k + i] += v_solution(i * h, k * tau);
           }
        }
    }
    /* if res(t=0)=0 the res(t=T) will be g_obs with sush rightside */
    SmoluchowskiSolver::forward_problem(smol_op, T, M, 
            RightSideGeneral(N, M, rightside), res);


    for(i = 0; i <= N; ++i) {
        /* g_obs = res(t=T) */
        ptmp[i] = res[(N + 1) * M + i];
        for(k = 0; k <= M; ++k) {
            /* res <- v0 */
            if((k < m1) || (m2 < k)) {
                res[(N + 1) * k + i] = 0;
            } else {
                res[(N + 1) * k + i] = v0(i * h, k * tau);
            }
        }
    }
    SourceFunctionGeneral v_res(N, H, M, T, m1, m2, res);
    StopCriteria stop_criteria(num_iter);
    CallbackPrintIter callback;
    solver.inverse_problem(v_res, ptmp, alpha, dzeta, stop_criteria, &callback); 



#ifdef CALC_DIRECT
    f = fopen("times_direct", "a");
#else
    f = fopen("times_fast", "a");
#endif
    fprintf(f, "N: %d, H: %f, M: %d, T: %f\n", N, H, M, T);
    fprintf(f, "apply linear operator max time: %f\n", smol_op.GetMaxApplyTime());
    fprintf(f, "apply adjoint linear operator max time: %f\n",
            smol_op.GetMaxApplyAdjointTime());
    fclose(f);

    delete [] ptmp;
    delete [] res;
    delete [] rightside;
    return 0;
}
