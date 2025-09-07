#include <smoluchowski/matrix.hpp>
#include <smoluchowski/cross.hpp>
#include <smoluchowski/solver.hpp>
#include <math.h>
#include <stdio.h>

double K(unsigned i, unsigned j, double h)
{
#if 1
    return pow(pow((i + 1) * h, 1.0 / 3.0) + pow((j + 1) * h, 1.0 / 3.0), 2) *
        sqrt(1.0 / ((i + 1) * h) + 1.0 / ((j + 1) * h));
#endif
    return 1.0;
    return pow((i + 1) * h * (j + 1) * h, 1.0 / 6.0);
#if 0
    return (pow((i + 1) * h, 1.0 / 3.0) + pow((j + 1) * h, 1.0 / 3.0)) *
        (1.0 / (pow((i + 1) * h, 1.0 / 3.0) * (j + 1) * h) +
         1.0 / (pow((j + 1) * h, 1.0 / 3.0) * (i + 1) * h));
#endif

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

int main(int argc, const char *const *argv)
{
    unsigned i, k;
    unsigned N, M;
    unsigned m1, m2;
    double H, T;
    double h;
    double *v0, *v_res, *v_solution, *c_obs_s, *c_obs_r;
    double cn, tmp;

    FILE *graph;
    graph = fopen("graph.csv", "r");
    if(!graph) {
        perror("graph.csv");
        fprintf(stderr, "start ./smolsolver_obs only after ./smolsolver\n");
        return 1;
    }
    fscanf(graph, "%u,%u,%lf,%lf,%u,%u\n", &N, &M, &H, &T, &m1, &m2);
    h = H / double(N);
    printf("N: %d, H: %f, M: %d, T: %f\n", N, H, M, T);
    printf("m1=%d, m2=%d\n", m1, m2);
    c_obs_s = new double[(N + 1) * (M + 1)];
    c_obs_r = new double[(N + 1) * (M + 1)];
    v0 = new double[(N + 1) * (M + 1)];
    v_res = new double[(N + 1) * (M + 1)];
    v_solution = new double[(N + 1) * (M + 1)];

    for(k = m1; k <= m2; ++k) {
        for(i = 0; i <= N; ++i) {
            tmp = fscanf(graph, "%lf,%lf,%lf\n", &v0[(N + 1) * k + i],
                    &v_res[(N + 1) * k + i], &v_solution[(N + 1) * k + i]);
            if(tmp != 3) {
                fprintf(stderr, "expected real value v0,v_res,v_sol\n");
                return 2;
            }
        }
    }
    fclose(graph);

    FunctionMatrix K_M(N + 1, N + 1, h, K);
    Cross K_approx(K_M);
    K_approx.Approximate(1e-6);
    FunctionMatrix Psi_M(N + 1, N + 1, h, Psi);
    Cross Psi_approx(Psi_M);
    Psi_approx.Approximate(1e-6);
    printf("K   rank: %d\n", K_approx.GetRank());
    printf("Psi rank: %d\n", Psi_approx.GetRank());

    SmoluchowskiCalcFast smol_base(K_approx, Psi_approx, H, N);
    SmoluchowskiNonLinearOperator smol_op_nl(smol_base);
    for(i = 0; i <= N; ++i) {
        c_obs_s[i] = c0(i * h);
        c_obs_r[i] = c0(i * h);
    }
    SmoluchowskiSolver::forward_problem(smol_op_nl, T, M,
            RightSideGeneral(N, M, v_solution), c_obs_s); 
    SmoluchowskiSolver::forward_problem(smol_op_nl, T, M, 
            RightSideGeneral(N, M, v_res), c_obs_r);

    FILE *graph_obs;
    graph_obs = fopen("graph_obs_new.csv", "w");
    if(!graph_obs) {
        perror("graph_obs_new.csv");
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
#if 0
    

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
    SmoluchowskiSolver::forward_problem(smol_op_nl, T * 50.0, M,
            RightSideGeneral(N, M, rightside), c_obs_s); 
    SmoluchowskiSolver::forward_problem(smol_op_nl, T * 50.0, M, 
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
    delete [] v0;
    delete [] v_res;
    delete [] v_solution;
    return 0;
}
