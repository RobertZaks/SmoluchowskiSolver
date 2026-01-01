#include <smoluchowski/matrix.hpp>
#include <smoluchowski/cross.hpp>
#include <smoluchowski/solver.hpp>
#include <math.h>
#include <stdio.h>

double K(unsigned i, unsigned j, double h)
{
    return pow((i + 1) * h * (j + 1) * h, 1.0);
}

double Psi(unsigned i, unsigned j, double h)
{
    return 0.0;
}

double c0(unsigned i, double h)
{
    if (i == 0)
        return exp(-h / 2.0) / (h / 2.0);
    else {
        double x = i * h;
        return exp(-x) / x;
    }
}

int main()
{
    unsigned i, k;
    unsigned N = 64000, M = 500;
    double H = 50.0, T = 2.0;
    double tau, h, x;
    double *res, *ptmp, *mass;
    h = H / double(N);
    tau = T / double(M);
    res = new double[(N + 1) * (M + 1)];
    ptmp = new double[N + 1];
    mass = new double[M + 1];
    
    
    printf("N: %d, H: %f, M: %d, T: %f\n", N, H, M, T);

    FunctionMatrix K_M(N + 1, N + 1, h, K);
    Cross K_approx(K_M);
    K_approx.Approximate(1e-6);
    FunctionMatrix Psi_M(N + 1, N + 1, h, Psi);
    Cross Psi_approx(Psi_M);
    Psi_approx.Approximate(1e-6);
    printf("K   rank: %d\n", K_approx.GetRank());
    printf("Psi rank: %d\n", Psi_approx.GetRank());

    for(i = 0; i <= N; ++i) {
        ptmp[i] = 0;
        res[i] = c0(i, h);
    }
    
    SmoluchowskiCalcFast smol_base(K_approx, Psi_approx, H, N);
    SmoluchowskiNonLinearOperator smol_op_nl(smol_base);
    SmoluchowskiSolver::forward_problem(smol_op_nl, T, M,
            RightSideNoTime(N, ptmp), res);

    FILE *graph;

    graph = fopen("smolgel_conc.csv", "w");
    for(i = 0; i <= N; i++) {
        x = double(i) * h;
        fprintf(graph, "%f\t%f\n", x, x * res[(N + 1) * M + i]);
    }
    fclose(graph);

    graph = fopen("smolgel_mass.csv", "w");
    for(k = 0; k <= M; k++) {
        mass[k] = 0;
        for(i = 0; i <= N; i++) {
            mass[k] += double(i) * res[(N + 1) * k + i];
        }
        mass[k] = h * h * (mass[k] - 0.5 * (0 + N * res[(N + 1) * k + N]));
        fprintf(graph, "%f\t%f\n", double(k) * tau, mass[k]);
    }
    fclose(graph);

    delete [] mass;
    delete [] ptmp;
    delete [] res;
    return 0;
}
