#include <smoluchowski/matrix.hpp>
#include <smoluchowski/cross.hpp>
#include <smoluchowski/solver.hpp>
#include <math.h>
#include <stdio.h>

/* This file contains the numerical method for solving
 * unsteady population balance equation with particles coagulation.
 * The calculated solution compared with analytical solution from
 * "An exact analytical solution to unsteady population balance
 * equation with particles coagulation" by
 * Eugenya V. Makoveeva, Dmitri V. Alexandrov,
 * https://doi.org/10.1016/j.cnsns.2024.107879.
 *
 * The max norm of difference on grid is printed.
 * To graph the results (if you are in the examples dir) you can run this
 * program by ``./partial`` and then run ``./scripts/draw_partial.py`` (you
 * need python and matplotlib to it; or you can draw it manually by gnuplot),
 * then open ``graph_partial.pdf``, ``graph_partial_m.pdf``.
 */

// coagulation kernel
double K(unsigned i, unsigned j, double h)
{
    return 1.0;
    return pow((i + 1) / double(j + 1), 1.0 / 3.0) +
        pow((j + 1) / double(i + 1), 1.0 / 3.0) + 2;
}

// fragmentation kernel
double Psi(unsigned i, unsigned j, double h)
{
    return 0.0;
}

// Calculate V(tau) - mass 
double calc_mass(unsigned N, double h, const double *g)
{
    unsigned i;
    double mass = 0;
    for(i = 0; i <= N; ++i) {
        mass += double(i) * g[i];
    }
    mass = h * h * (mass - 0.5 * double(N) * g[N]);
    return mass;
}

double calc_density(unsigned N, double h, const double *g)
{
    unsigned i;
    double dens = 0;
    for(i = 0; i <= N; ++i) {
        dens += g[i];
    }
    dens = h * (dens - 0.5 * (g[0] + g[N]));
    return dens;
}

int main()
{
    /* grid */
    // we assume supp(K)\subset[0;H]^{2}
    double H = 40, T = 20.0;
    unsigned N = 4000, M = unsigned(T * pow(H / double(N), -2));
    double h = H / double(N), tau = T / double(M);
 
    /* physical constants */
    double cs = 10, kappa = 0.2, chi = 1e-2, gamma = 1.0;
    double phi00 = 1.0, b0 = 1, delta0 = 0.2;

    double q_ci_cs = delta0 + 1 + phi00 / (cs * b0 * b0);
    
    /* Absorbing layers */
    // Ns - number of absorbing layers
    unsigned Ns = 0;
    // the decreasing rate
    double d = 5.0;

    unsigned Nt = N + Ns;

    printf("N: %d, H: %f, Ns: %d, M: %d, T: %f\n", N, H, Ns, M, T);

    if(N < 5) {
        fputs("Too small nodes on grid for calc; exiting now\n", stderr);
        return 1;
    }
    
    /* approximate kernels */
    double tol_cross = 1e-6;
    FunctionMatrix K_M(N + 1, N + 1, h, K);
    Cross K_approx(K_M);
    K_approx.Approximate(tol_cross);
    printf("K   rank: %d\n", K_approx.GetRank());
    FunctionMatrix Psi_M(N + 1, N + 1, h, Psi);
    Cross Psi_approx(Psi_M);
    Psi_approx.Approximate(tol_cross);
    /* init smoluchowski calculation base class */
    SmoluchowskiCalcFast smol_base(K_approx, Psi_approx, H, N);
    /* create nonlinear smoluchowski operator:
     * -0.5 * L1[g,g] + g * L5[g] */
    SmoluchowskiNonLinearOperator smol_op_nl(smol_base);
    
    double delta_g, h2 = pow(h, 2.0);
    double init_mass = phi00/(b0 * b0);
    double *phi_old, *phi_new, *ptmp;
    double *mass = new double[M + 1];
    double *density = new double[M + 1];

    phi_old = new double[Nt + 1];
    phi_new = new double[Nt + 1];
    double tmp = 0.0;
#if 1
    for(unsigned i = 0; i <= Nt; ++i) {
        phi_new[i] = phi00 * exp(-b0 * double(i) * h);
    }
#else
    double mean = 5.0, std = 2.0;
    for(unsigned i = 0; i <= Nt; ++i) {
        //phi_new[i] = exp(-pow((double(i) * h - mean) / std, 2));
        phi_new[i] = phi00 * ( 2./sqrt(M_PI)/b0 *exp(-double(i*i) * h * h));
    }
#endif
    tmp = calc_mass(Nt, h, phi_new);
    for(unsigned i = 0; i <= Nt; ++i) {
        phi_new[i] *= init_mass / tmp;
    }
    printf("mass_cmp:%f\n", init_mass - calc_mass(Nt, h, phi_new));
    /* calc phi */
    for(unsigned k = 0; k <= M; ++k) {
        ptmp = phi_new;
        phi_new = phi_old;
        phi_old = ptmp;
        smol_op_nl.Apply(phi_old, phi_new);
        mass[k] = calc_mass(Nt, h, phi_old);
        density[k] = calc_density(Nt, h, phi_old);
        delta_g = pow(q_ci_cs - 1 - mass[k] / cs, gamma);
        if(k == M) {
            printf("delta_g=%f\n", delta_g);
        }
        phi_new[0] += delta_g * (
            kappa * (
                -1.5 * phi_old[0] +
                2 * phi_old[1] +
                -0.5 * phi_old[2]
            ) / h -
            chi * (
                2 * phi_old[0] +
                -5 * phi_old[1] +
                4 * phi_old[2] +
                -phi_old[3]
            ) / h2
        );
        for(unsigned i = 1; i < Nt; ++i) {
            phi_new[i] += delta_g * (
                kappa * (phi_old[i + 1] - phi_old[i - 1]) * 0.5 / h -
                chi * (phi_old[i + 1] - 2 * phi_old[i] + phi_old[i - 1]) / h2
            );
        }
        phi_new[Nt] += delta_g * (
            kappa * (
                0.5 * phi_old[Nt - 2] +
                -2 * phi_old[Nt - 1] +
                1.5 * phi_old[Nt]
            ) / h -
            chi * (
                -phi_old[Nt - 3] +
                4 * phi_old[Nt - 2] +
                -5 * phi_old[Nt - 1] +
                2 * phi_old[Nt]
            ) / h2
        );
        for(unsigned i = 0; i <= Nt; ++i) {
            phi_new[i] = phi_old[i] - tau * phi_new[i];
        }
        for(unsigned i = N + 1; i <= Nt; ++i) {
            phi_new[i] *= exp(-d * double(i - N) * h);
        }
    }

    unsigned Nb = 40000;
    double b, b1 = 0;
    double b_step = (b1 - b0) / double(Nb);
    double *delta = new double[Nb + 1];
    for(unsigned k = 0; k <= Nb; ++k) {
        b = b0 + b_step * k;
        if(gamma == 1.0) {
            delta[k] = delta0 * pow(b / b0, 2 * kappa / cs) * 
                exp(2 * chi * (b - b0) / cs);
        } else {
            delta[k] = pow(pow(delta0, 1 - gamma) + kappa * (1 - gamma) / cs *
            (log(pow(b, 2) / pow(b0, 2)) + 2 * chi * (b - b0) / kappa),
            1.0 / (1.0 - gamma));
        }
    }
    double *hb = new double[Nb + 1];
    double hb_int = 0.5 * pow(delta[0], gamma) * (kappa + chi * b0) / b0;
    for(unsigned k = 0; k <= Nb; ++k) {
        b = b0 + b_step * k;
        tmp = pow(delta[k], gamma) * (kappa + chi * b) / b;
        hb_int += tmp;
        hb[k] = -0.5 * phi00 * pow(b / b0, 2) +
            b * b * b_step * (hb_int - 0.5 * tmp);
    }

    unsigned M_tau = 0;
    double *taub = new double[Nb + 1];
    hb_int = 0.5 / hb[0];
    for(unsigned k = 0; k <= Nb; ++k) {
        tmp = 1.0 / hb[k];
        hb_int += tmp;
        taub[k] = b_step * (hb_int - 0.5 * tmp);
        if(taub[k] >= T) {
            M_tau = k;
            printf("tau=%f\n", taub[k]);
            break;
        }
    }

    FILE *graph;
    double *analytic_solution = new double[N + 1];
    graph = fopen("graph_partial_analytic.csv", "w");
    fprintf(graph, "%f\t%d\t%d\t%f\t%d\n", H, N, Ns, T, M);
    printf("delta[tau]=%f\n", delta[M_tau]);
    for(unsigned i = 0; i <= N; ++i) {
        b = b0 + b_step * M_tau;
        analytic_solution[i] =
                -2 * hb[M_tau] * exp(-b * h * double(i));
        fprintf(graph, "%f\t%f\n", double(i) * h,
                analytic_solution[i]);
    }
    fclose(graph);

    graph = fopen("graph_partial_solution.csv", "w");
    fprintf(graph, "%f\t%d\t%d\t%f\t%d\n", H, N, Ns, T, M);
    for(unsigned i = 0; i <= Nt; ++i) {
        fprintf(graph, "%f\t%f\n", double(i) * h, phi_new[i]);
    }
    fclose(graph);

    graph = fopen("graph_partial_solution_mass.csv", "w");
    //fprintf(graph, "%f\t%d\t%d\t%f\t%d\n", H, N, Ns, T, M);
    for(unsigned k = 0; k <= M; ++k) {
        fprintf(graph, "%E\t %E\t %E\n",
                double(k) * tau, mass[k], density[k]);
    }
    fclose(graph);

    double cn = -1;
    for(unsigned i = 0; i <= N; ++i) {
        tmp = fabs(phi_new[i] - analytic_solution[i]);
        if(cn < tmp) {
            cn = tmp;
        }
    }
    printf("C-norm:%e\n", cn);

    delete [] analytic_solution;
    delete [] taub;
    delete [] delta;
    delete [] hb;
    delete [] phi_old;
    delete [] phi_new;
    return 0;
}
