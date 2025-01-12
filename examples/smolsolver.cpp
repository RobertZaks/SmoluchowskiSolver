#include <smoluchowski/matrix.hpp>
#include <smoluchowski/cross.hpp>
#include <smoluchowski/solver.hpp>
#include <math.h>
#include <stdio.h>

double K(unsigned i, unsigned j, double h)
{
#if 0
    return (pow((i + 1) * h, 1.0 / 3.0) + pow((j + 1) * h, 1.0 / 3.0)) *
        (1.0 / (pow((i + 1) * h, 1.0 / 3.0) * (j + 1) * h) +
         1.0 / (pow((j + 1) * h, 1.0 / 3.0) * (i + 1) * h));
#endif
#if 0
    return pow(pow((i + 1) * h, 1.0 / 3.0) + pow((j + 1) * h, 1.0 / 3.0), 2) *
        sqrt(1.0 / ((i + 1) * h) + 1.0 / ((j + 1) * h));
#endif
    return 1.0;
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
    //return exp(-x);
    //return exp(-2.0 * x);
#if 0
    return exp(-(x - 1.5) * (x - 1.5) / 2.0) +
        2.0 * exp(-(x - 4.5) * (x - 4.5) / 2.0);
#endif
    return (exp(-(x - 1.5) * (x - 1.5) / 2.0) +
        2.0 * exp(-(x - 4.5) * (x - 4.5) / 2.0)) * 
        tanh(4.0 * t);
}

double v0(double x, double t)
{
    //return exp(-2.0 * x);
    //return exp(-x);
    return exp(-(x - 3.0) * (x - 3.0) / 4.0);
}

class CallbackPrintIterSaveFunc : public SmoluchowskiSolver::Callback {
    FILE *vals_func;
public:
    CallbackPrintIterSaveFunc(const char *name) {
        vals_func = fopen(name, "w");
    }
    ~CallbackPrintIterSaveFunc() {
        fclose(vals_func);
    }
    void InverseProblemIter(const StopCriteria& stop_criteria) {
        fprintf(stderr, "iter: %d\n", stop_criteria.GetCurIter());
        if(stop_criteria.NeedCalcFunc()) {
            fprintf(vals_func, "%e\n", stop_criteria.GetCurFunc());
        }
    }
};

int main()
{
    unsigned i, k;
    unsigned N = 80, M = 80;
    unsigned m1, m2;
    double H = 6.0, T = 0.5;
    double t1, t2, tau, h, alpha, dzeta, c1, tol = 1e-6;
    double *res, *rightside, *ptmp, *c_obs_s, *c_obs_r;
    double cn, tmp;
    h = H / double(N);
    tau = T / double(M);
    rightside = new double[(N + 1) * (M + 1)];
    res = new double[(N + 1) * (M + 1)];
    ptmp = new double[N + 1];
    c_obs_s = new double[(N + 1) * (M + 1)];
    c_obs_r = new double[(N + 1) * (M + 1)];
    
    c1 = T;
    alpha = 0.0001;
    dzeta = 2.0 / (2.0 * alpha + c1);
    t1 = T / 5;
    t2 = 3 * T / 5;
    m1 = floor(t1 / T * M);
    m2 = ceil(t2 / T * M);
    printf("N: %d, H: %f, M: %d, T: %f\n", N, H, M, T);
    printf("m1=%d, m2=%d\n", m1, m2);
    printf("c1=%f\n", c1);

    FunctionMatrix K_M(N + 1, N + 1, h, K);
    Cross K_approx(K_M);
    K_approx.Approximate(1e-6);
    FunctionMatrix Psi_M(N + 1, N + 1, h, Psi);
    Cross Psi_approx(Psi_M);
    Psi_approx.Approximate(1e-6);
    printf("K   rank: %d\n", K_approx.GetRank());
    printf("Psi rank: %d\n", Psi_approx.GetRank());

    for(i = 0; i <= N; ++i) {
        ptmp[i] = c0(i * h);
        res[i] = 0;
    }
    
    SmoluchowskiCalcFast smol_base(K_approx, Psi_approx, H, N);
    SmoluchowskiLinearOperator smol_op(smol_base, ptmp);
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
    unsigned num_iter = 7;
    //StopCriteria stop_criteria(num_iter, true);
    StopCriteria stop_criteria(tol);
    CallbackPrintIterSaveFunc callback("functional.csv");
    solver.inverse_problem(v_res, ptmp, alpha, dzeta, stop_criteria, &callback);

    FILE *graph;
    graph = fopen("graph.csv", "w");
    if(!graph) {
        perror("graph.csv");
    }
    fprintf(graph, "%d,%d,%f,%f,%f,%f\n", N, M, H, T, t1, t2);
    cn = 0;
    for(k = m1; k <= m2; ++k) {
        for(i = 0; i <= N; ++i) {
            fprintf(graph, "%f,%f,%f\n", v0(i * h, k * tau), 
                    v_res(i, k), v_solution(i * h, k * tau));
            tmp = fabs(v_res(i, k) - v_solution(i * h, k * tau));
            if(cn < tmp) {
                cn = tmp;
            }
        }
    }
    fclose(graph);
    printf("cnorm=%e\n", cn);


    SmoluchowskiNonLinearOperator smol_op_nl(smol_base);
    for(i = 0; i <= N; ++i) {
        for(k = 0; k <= M; ++k) {
            if((k < m1) || (m2 < k)) {
                rightside[(N + 1) * k + i] = 0;
            } else {
                rightside[(N + 1) * k + i] = v_solution(i * h, k * tau);
            }
        }
    }
    for(i = 0; i <= N; ++i) {
        c_obs_s[i] = c0(i * h);
        c_obs_r[i] = c0(i * h);
    }
    SmoluchowskiSolver::forward_problem(smol_op_nl, T, M,
            RightSideGeneral(N, M, rightside), c_obs_s); 
    SmoluchowskiSolver::forward_problem(smol_op_nl, T, M, 
            RightSideSourceFunction(v_res), c_obs_r);

    FILE *graph_obs;
    graph_obs = fopen("graph_obs.csv", "w");
    if(!graph) {
        perror("graph_obs.csv");
    }
    cn = -1;
    for(i = 0; i <= N; ++i) {
        fprintf(graph_obs, "%f,%f\n", 
                c_obs_s[(N + 1) * M + i], c_obs_r[(N + 1) * M + i]);
        tmp = fabs(c_obs_r[(N + 1) * M + i] - c_obs_s[(N + 1) * M + i]);
        if(cn < tmp) {
            cn = tmp;
        }
    }
    fclose(graph_obs);
    printf("cnorm_obs=%e\n", cn);

    /* fixtime: */
#if 1
    

    for(i = 0; i <= N; ++i) {
        res[i] = v0(i * h, m1);
    }
    for(i = 0; i <= M; ++i) {
        res[(N + 1) + i] = 1;
    }
    SourceFunctionFixTime v_res_fixtime(N, H, M, T, m1, m2, 
            res + (N + 1), res);
    stop_criteria.Reset();
    CallbackPrintIterSaveFunc callback_fixtime("functional_fixtime.csv");
    solver.inverse_problem(v_res_fixtime, ptmp, alpha, dzeta, stop_criteria,
            &callback_fixtime);

    graph = fopen("graph_fixtime.csv", "w");
    if(!graph) {
        perror("graph_fixtime.csv");
    }
    cn = -1;
    for(i = 0; i <= N; ++i) {
        fprintf(graph, "%f,%f,%f\n", v0(i * h, m1), 
            v_res_fixtime(i, m1), v_solution(i * h, m1));
        tmp = fabs(v_res_fixtime(i, m1) - v_solution(i * h, m1));
        if(cn < tmp) {
            cn = tmp;
        }
    }
    fclose(graph);
    printf("cnorm fix_time=%e\n", cn);

    for(i = 0; i <= N; ++i) {
        for(k = 0; k <= M; ++k) {
            if((k < m1) || (m2 < k)) {
                rightside[(N + 1) * k + i] = 0;
            } else {
                rightside[(N + 1) * k + i] = v_solution(i * h, k * tau);
            }
        }
    }
    for(i = 0; i <= N; ++i) {
        c_obs_s[i] = c0(i * h);
        c_obs_r[i] = c0(i * h);
    }
    SmoluchowskiSolver::forward_problem(smol_op_nl, T, M,
            RightSideGeneral(N, M, rightside), c_obs_s); 
    SmoluchowskiSolver::forward_problem(smol_op_nl, T, M, 
            RightSideSourceFunction(v_res_fixtime), c_obs_r);

    graph_obs = fopen("graph_obs_fixtime.csv", "w");
    if(!graph) {
        perror("graph_obs_fixtime.csv");
    }
    cn = -1;
    for(i = 0; i <= N; ++i) {
        fprintf(graph_obs, "%f,%f\n",
                c_obs_s[(N + 1) * M + i], c_obs_r[(N + 1) * M + i]);
        tmp = fabs(c_obs_r[(N + 1) * M + i] - c_obs_s[(N + 1) * M + i]);
        if(cn < tmp) {
            cn = tmp;
        }
    }
    fclose(graph_obs);
    printf("cnorm_obs_fixtime=%e\n", cn);
#endif

    delete [] c_obs_r;
    delete [] c_obs_s;
    delete [] ptmp;
    delete [] res;
    delete [] rightside;
    return 0;
}
