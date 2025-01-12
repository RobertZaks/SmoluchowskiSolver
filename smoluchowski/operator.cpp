#include "operator.hpp"
#include <math.h>

SmoluchowskiCalc::SmoluchowskiCalc(const Skeleton& K_skel,
    const Skeleton& Psi_skel, double _H, unsigned _N)
    : K(K_skel), Psi(Psi_skel), N(_N), H(_H), h(H / double(N))
{ }


/*=============SmoluchowskiCalcDirect===============*/

SmoluchowskiCalcDirect::SmoluchowskiCalcDirect(const Skeleton& K_skel,
    const Skeleton& Psi_skel, double _H, unsigned _N)
    : SmoluchowskiCalc(K_skel, Psi_skel, _H, _N), Uf(0)
{ }

SmoluchowskiCalcDirect::~SmoluchowskiCalcDirect()
{
    if(Uf) {
        delete [] Uf;
    }
}

//TODO: no use of GetU and GetV it direct calc!
void SmoluchowskiCalcDirect::calc_L1(const double *f, 
        const double *g, double *res)
{
    unsigned i, j, k, r;
    r = K.GetRank();
    for(i = 0; i <= N; ++i) {
        res[i] = 0;
        for(k = 0; k < r; ++k) {
            for(j = 0; j <= i; ++j) {
                res[i] += K.GetU().elem(i - j, k) * f[i - j] *
                    K.GetV().elem(j, k) * g[j];
            }
            res[i] -= 0.5 * (K.GetU().elem(i, k) * f[i] *
                    K.GetV().elem(0, k) * g[0] +
                    K.GetU().elem(0, k) * f[0] *
                    K.GetV().elem(i, k) * g[i]);
        }
        res[i] *= h;
    }
}

void SmoluchowskiCalcDirect::calc_L2(const double *f, 
        const double *g, double *res)
{
    double tmp;
    unsigned i, j, k, r;
    const Matrix& U = K.GetU();
    const Matrix& V = K.GetV();
    r = K.GetRank();
    for(i = 0; i <= N; ++i) {
        res[i] = 0;
        for(k = 0; k < r; ++k) {
            tmp = 0;
            for(j = i; j <= N; ++j) {
                tmp += U.elem(j - i, k) * f[j - i] * g[j];
            }
            tmp -= 0.5 * (U.elem(0, k) * f[0] * g[i] + 
                    U.elem(N - i, k) * f[N - i] * g[N]);
            res[i] += V.elem(i, k) * tmp;
        }
        res[i] *= h;
    }
}

void SmoluchowskiCalcDirect::fix_f_arg(const double *f)
{
    unsigned i, k, r;
    r = K.GetRank();
    if(!Uf) {
        Uf = new double[(N + 1) * K.GetRank()];
    }
    for(k = 0; k < r; ++k) {
        for(i = 0; i <= N; ++i) {
            Uf[(N + 1) * k + i] = K.GetU().elem(i, k) * f[i];
        }
    }
}

void
SmoluchowskiCalcDirect::calc_L1_fix_f(const double *g, double *res)
{
    unsigned i, j, k, r;
#if HANDLEXCEPTION
    if(!Uf) {
        throw "f must be set by fix_f_arg";
    }
#endif
    r = K.GetRank();
    for(i = 0; i <= N; ++i) {
        res[i] = 0;
        for(k = 0; k < r; ++k) {
            for(j = 0; j <= i; ++j) {
                res[i] += Uf[(N + 1) * k + i - j] * 
                    K.GetV().elem(j, k) * g[j];
            }
            res[i] -= 0.5 * (Uf[(N + 1) * k + i] * 
                    K.GetV().elem(0, k) * g[0] +
                    Uf[(N + 1) * k] * 
                    K.GetV().elem(i, k) * g[i]);
        }
        res[i] *= h;
    }
}

void
SmoluchowskiCalcDirect::calc_L2_fix_f(const double *g, double *res)
{
    double tmp;
    unsigned i, j, k, r;
#if HANDLEXCEPTION
    if(!Uf) {
        throw "f must be set by fix_f_arg";
    }
#endif
    r = K.GetRank();
    for(i = 0; i <= N; ++i) {
        res[i] = 0;
        for(k = 0; k < r; ++k) {
            tmp = 0;
            for(j = i; j <= N; ++j) {
                tmp += Uf[(N + 1) * k + j - i] * g[j];
            }
            tmp -= 0.5 * (Uf[(N + 1) * k] * g[i] +
                    Uf[(N + 1) * k + N - i] * g[N]);
            res[i] += K.GetV().elem(i, k) * tmp;
        }
        res[i] *= h;
    }
}



void SmoluchowskiCalcDirect::calc_L3(const double *g, double *res) const
{
    unsigned i, j;
    for(i = 0; i <= N; ++i) {
        res[i] = 0;
        for(j = 0; j <= i; ++j) {
            res[i] += Psi.GetM().elem(i, j) * g[j];
        }
        res[i] = h * (res[i] - 0.5 * (Psi.GetM().elem(i, 0) * g[0] +
                    Psi.GetM().elem(i, i) * g[i]));
    }
}

void SmoluchowskiCalcDirect::calc_L4(const double *g, double *res) const
{
    unsigned i, j;
    for(i = 0; i <= N; ++i) {
        res[i] = 0;
        for(j = i; j <= N; ++j) {
            res[i] += Psi.GetM().elem(j, i) * g[j];
        }
        res[i] = h * (res[i] - 0.5 * (Psi.GetM().elem(i, i) * g[i] + 
                    Psi.GetM().elem(N, i) * g[N]));
    }
}

void SmoluchowskiCalcDirect::calc_L5(const double *g, double *res) const
{
    unsigned i, j;
    for(i = 0; i <= N; ++i) {
        res[i] = 0;
        for(j = 0; j <= N; ++j) {
            res[i] += K.GetM().elem(i, j) * g[j];
        }
        res[i] = h * (res[i] - 0.5 * (K.GetM().elem(i, 0) * g[0] + 
                    K.GetM().elem(i, N) * g[N]));
    }
}

/*=============SmoluchowskiCalcFast===============*/


SmoluchowskiCalcFast::SmoluchowskiCalcFast(const Skeleton& K_skel,
    const Skeleton& Psi_skel, double _H, unsigned _N)
    : SmoluchowskiCalc(K_skel, Psi_skel, _H, _N), 
    N_log2(ceil(log2(2 * N + 1))), N_2(1 << N_log2), Fourier(N_log2),
    Uf(0), Uf_fourier(0)
{
    unsigned r1, r2;
    r1 = K.GetRank();
    r2 = Psi.GetRank();
    tmp_mv = new double[r1 > r2 ? r1 : r2];
    tmp_complex0 = new Complex[N_2];
    tmp_complex1 = new Complex[N_2];
    tmp_complex2 = new Complex[N_2];
}

SmoluchowskiCalcFast::~SmoluchowskiCalcFast()
{
    delete [] tmp_mv;
    delete [] tmp_complex0;
    delete [] tmp_complex1;
    delete [] tmp_complex2;
    if(Uf) {
        delete [] Uf;
        delete [] Uf_fourier;
    }
}


void SmoluchowskiCalcFast::calc_L1(const double *f, 
        const double *g, double *res)
{
    unsigned i, k, r;
    r = K.GetRank();
    for(i = 0; i <= N; ++i) {
        res[i] = 0;
    }
    for(k = 0; k < r; ++k) {
        for(i = 0; i <= N; ++i) {
            tmp_complex0[i] = K.GetU().elem(i, k) * f[i];
        }
        Fourier.MatvecFast(tmp_complex0, tmp_complex1);
        for(i = 0; i <= N; ++i) {
            tmp_complex0[i] = K.GetV().elem(i, k) * g[i];
        }
        Fourier.MatvecFast(tmp_complex0, tmp_complex2);
        for(i = 0; i < N_2; ++i) {
            tmp_complex1[i] *= tmp_complex2[i];
        }
        Fourier.MatvecFast(tmp_complex1, tmp_complex2, 
                FourierTransform::conj);
        for(i = 0; i <= N; ++i) {
            res[i] +=
                h * (tmp_complex2[i].GetRe() / double(N_2) - 
                        0.5 * (K.GetU().elem(i, k) * f[i] *
                            tmp_complex0[0].GetRe() +
                            K.GetU().elem(0, k) * f[0] *
                            tmp_complex0[i].GetRe()));
        }
    }
}

void SmoluchowskiCalcFast::calc_L2(const double *f, 
        const double *g, double *res)
{
    unsigned i, k, r;
    r = K.GetRank();
    for(i = 0; i <= N; ++i) {
        res[i] = 0;
    }
    for(k = 0; k < r; ++k) {
        for(i = 0; i <= N; ++i) {
            tmp_complex0[N_2 - 1 - i] = g[i];
        }
        Fourier.MatvecFast(tmp_complex0, tmp_complex1);
        for(i = 0; i <= N; ++i) {
            tmp_complex0[i] = K.GetU().elem(i, k) * f[i];
        }
        for(i = N + 1; i < N_2; ++i) {
            tmp_complex0[i] = 0;
        }
        Fourier.MatvecFast(tmp_complex0, tmp_complex2);
        for(i = 0; i < N_2; ++i) {
            tmp_complex1[i] *= tmp_complex2[i];
        }
        Fourier.MatvecFast(tmp_complex1, tmp_complex2,
                FourierTransform::conj);
        for(i = 0; i <= N; ++i) {
            res[i] += h * K.GetV().elem(i, k) *
                (tmp_complex2[N_2 - 1 - i].GetRe() / double(N_2) -
                 0.5 * (tmp_complex0[0].GetRe() * g[i] +
                     tmp_complex0[N - i].GetRe() * g[N]));
        }
    }
}

void SmoluchowskiCalcFast::fix_f_arg(const double *f)
{
    unsigned i, k, r;
    r = K.GetRank();
    if(!Uf) {
        Uf = new double[(N + 1) * K.GetRank()];
        Uf_fourier = new Complex[N_2 * r];
    }
    for(k = 0; k < r; ++k) {
        for(i = 0; i <= N; ++i) {
            Uf[(N + 1) * k + i] = K.GetU().elem(i, k) * f[i];
            tmp_complex1[i] = Uf[(N + 1) * k + i];
        }
        Fourier.MatvecFast(tmp_complex1, Uf_fourier + N_2 * k);
    }
}

void
SmoluchowskiCalcFast::calc_L1_fix_f(const double *g, double *res)
{
    unsigned i, k, r;
#if HANDLEXCEPTION
    if(!Uf) {
        throw "f must be set by fix_f_arg";
    }
#endif
    r = K.GetRank();
    for(i = 0; i <= N; ++i) {
        res[i] = 0;
    }
    for(k = 0; k < r; ++k) {
        for(i = 0; i <= N; ++i) {
            tmp_complex0[i] = K.GetV().elem(i, k) * g[i];
        }
        Fourier.MatvecFast(tmp_complex0, tmp_complex1);
        for(i = 0; i < N_2; ++i) {
            tmp_complex1[i] *= Uf_fourier[N_2 * k + i];
        }
        Fourier.MatvecFast(tmp_complex1, tmp_complex2, 
                FourierTransform::conj);
        for(i = 0; i <= N; ++i) {
            res[i] += h * (tmp_complex2[i].GetRe() / double(N_2) - 
                    0.5 * (Uf[(N + 1) * k + i] * tmp_complex0[0].GetRe() +
                        Uf[(N + 1) * k] * tmp_complex0[i].GetRe()));
        }
    }
}

void
SmoluchowskiCalcFast::calc_L2_fix_f(const double *g, double *res)
{
    unsigned i, k, r;
#if HANDLEXCEPTION
    if(!Uf) {
        throw "f must be set by fix_f_arg";
    }
#endif
    r = K.GetRank();
    for(i = 0; i <= N; ++i) {
        res[i] = 0;
    }
    for(k = 0; k < r; ++k) {
        for(i = 0; i <= N; ++i) {
            tmp_complex0[N_2 - 1 - i] = g[i];
        }
        Fourier.MatvecFast(tmp_complex0, tmp_complex1);
        for(i = 0; i < N_2; ++i) {
            tmp_complex1[i] *= Uf_fourier[N_2 * k + i];
            tmp_complex0[i] = 0; // we must keep last part with zeros to K1
        }
        Fourier.MatvecFast(tmp_complex1, tmp_complex2,
                FourierTransform::conj);
        for(i = 0; i <= N; ++i) {
            res[i] += h * K.GetV().elem(i, k) *
                (tmp_complex2[N_2 - 1 - i].GetRe() / double(N_2) -
                 0.5 * (Uf[(N + 1) * k] * g[i] +
                     Uf[(N + 1) * k + N - i] * g[N]));
        }
    }
}


void SmoluchowskiCalcFast::calc_L3(const double *g, double *res) const
{
    unsigned i;
    Psi.MatvecLT(h, g, res);
    for(i = 0; i <= N; ++i) {
        res[i] -= 0.5 * h * 
            (Psi.GetM().elem(i, 0) * g[0] + Psi.GetM().elem(i, i) * g[i]);
    }
}

void SmoluchowskiCalcFast::calc_L4(const double *g, double *res) const
{
    unsigned i;
    Psi.MatvecUT(h, g, res, Skeleton::transpose);
    for(i = 0; i <= N; ++i) {
        res[i] -= 0.5 * h * 
            (Psi.GetM().elem(i, i) * g[i] + Psi.GetM().elem(N, i) * g[N]);
    }
}

void SmoluchowskiCalcFast::calc_L5(const double *g, double *res) const
{
    unsigned i;
    K.MatvecF(h, g, res, tmp_mv);
    for(i = 0; i <= N; ++i) {
        res[i] -= 0.5 * h * 
            (K.GetM().elem(i, 0) * g[0] + K.GetM().elem(i, N) * g[N]);
    }
}


/*=============SmoluchowskiOperator===============*/

SmoluchowskiOperator::SmoluchowskiOperator(SmoluchowskiCalc& calc)
    : smol_base(calc)
{
    N = smol_base.GetN();
    H = smol_base.GetH();
    h = H / double(N);
    tmp_l = new double[N + 1];
}

SmoluchowskiOperator::~SmoluchowskiOperator()
{
    delete [] tmp_l;
}

/*=============SmoluchowskiNonLinearOperator===============*/

SmoluchowskiNonLinearOperator::SmoluchowskiNonLinearOperator(
        SmoluchowskiCalc& calc)
    : SmoluchowskiOperator(calc)
{ }

SmoluchowskiNonLinearOperator::~SmoluchowskiNonLinearOperator()
{ }

void SmoluchowskiNonLinearOperator::Apply(const double *g, double *res)
{
    unsigned i;
    for(i = 0; i <= N; ++i) {
        res[i] = double(i) * h;
    }
    smol_base.calc_L3(res, tmp_l);
    smol_base.calc_L4(g, res);
    for(i = 1; i <= N; ++i) {
        res[i] -= g[i] / (double(i) * h) * tmp_l[i];
    }
    smol_base.calc_L5(g, tmp_l);
    for(i = 0; i <= N; ++i) {
        res[i] -= g[i] * tmp_l[i];
    }
    smol_base.calc_L1(g, g, tmp_l);
    for(i = 0; i <= N; ++i) {
        res[i] = -(res[i] + 0.5 * tmp_l[i]);
    }
}

/*=============SmoluchowskiLinearOperator===============*/

SmoluchowskiLinearOperator::SmoluchowskiLinearOperator(
        SmoluchowskiCalc& calc, const double *st)
    : SmoluchowskiOperator(calc)
{
    unsigned i;
    c0 = new double[N + 1];
    for(i = 0; i <= N; ++i) {
        c0[i] = st[i];
    }
    smol_base.fix_f_arg(c0);
    a = new double[N + 1];
    calc_a();
}

SmoluchowskiLinearOperator::~SmoluchowskiLinearOperator()
{
    delete [] c0;
    delete [] a;
}

const double *SmoluchowskiLinearOperator::Get_c0() const
{
    return c0;
}

const double *SmoluchowskiLinearOperator::Get_a() const
{
    return a;
}


void SmoluchowskiLinearOperator::calc_a()
{
    unsigned i;
    for(i = 0; i <= N; ++i) {
        a[i] = double(i) * h;
    }
    smol_base.calc_L3(a, tmp_l);
    smol_base.calc_L5(c0, a);
    for(i = 1; i <= N; ++i) {
        a[i] += tmp_l[i] / ((double)(i) * h);
    }
}

void SmoluchowskiLinearOperator::Apply(const double *g, double *res)
{
    unsigned i;
    smol_base.calc_L5(g, res);
    smol_base.calc_L1_fix_f(g, tmp_l);
    for(i = 0; i <= N; ++i) {
        res[i] = c0[i] * res[i] - tmp_l[i];
    }
    smol_base.calc_L4(g, tmp_l);
    for(i = 0; i <= N; ++i) {
        res[i] += a[i] * g[i] - tmp_l[i];
    }
}

void SmoluchowskiLinearOperator::ApplyAdjoint(const double *g, double *res)
{
    unsigned i;
    for(i = 0; i <= N; ++i) {
        tmp_l[i] = c0[i] * g[i];
    }
    smol_base.calc_L5(tmp_l, res);
    smol_base.calc_L2_fix_f(g, tmp_l);
    for(i = 0; i <= N; ++i) {
        res[i] -= tmp_l[i];
    }
    smol_base.calc_L3(g, tmp_l);
    for(i = 0; i <= N; ++i) {
        res[i] += a[i] * g[i] - tmp_l[i];
    }
}
