#include <smoluchowski/matrix.hpp>
#include <smoluchowski/cross.hpp>
#include <smoluchowski/solver.hpp>
#include <math.h>
#include <stdio.h>

double K(unsigned i, unsigned j, double h)
{
    return 1.0;
}

double Psi(unsigned i, unsigned j, double h)
{
    return 0.0;
}

double c0(double x)
{
    return exp(-x);
}

double analytic_solution(double x, double t)
{
    return exp(-x / (1 + 0.5 * t)) * pow(1 + 0.5 * t, -2);
}

int main()
{
    unsigned i;
    unsigned N = 640, M = 640;
    double H = 12.0, T = 3.0;
    double tau, h, x;
    double *res, *ptmp;
    double cn, tmp, tol_cross = 1e-6;
    h = H / double(N);
    tau = T / double(M);
    res = new double[(N + 1) * (M + 1)];
    ptmp = new double[N + 1];
    
    printf("N: %d, H: %f, M: %d, T: %f\n", N, H, M, T);

    FunctionMatrix K_M(N + 1, N + 1, h, K);
    Cross K_approx(K_M);
    K_approx.Approximate(tol_cross);
    FunctionMatrix Psi_M(N + 1, N + 1, h, Psi);
    Cross Psi_approx(Psi_M);
    Psi_approx.Approximate(tol_cross);

    for(i = 0; i <= N; ++i) {
        ptmp[i] = 0;
        res[i] = c0(i * h);
    }
    
    SmoluchowskiCalcFast smol_base(K_approx, Psi_approx, H, N);
    SmoluchowskiNonLinearOperator smol_op_nl(smol_base);
    SmoluchowskiSolver::forward_problem(smol_op_nl, T, M,
            RightSideNoTime(N, ptmp), res);

    FILE *graph;
    graph = fopen("graph_smoltest.csv", "w");
    cn = -1;
    for(i = 0; i <= N; ++i) {
        x = double(i) * h;
        fprintf(graph, "%f\t%f\t%f\n", x,
                x * res[(N + 1) * M + i],
                x * analytic_solution(double(i) * h, double(M) * tau));
        tmp = fabs(res[(N + 1) * M + i] -
                analytic_solution(double(i) * h, double(M) * tau));
        if(cn < tmp) {
            cn = tmp;
        }
    }
    printf("C-norm:%e\n", cn);
    fclose(graph);

    delete [] ptmp;
    delete [] res;
    return 0;
}
