#include "skeleton.hpp"

Skeleton::Skeleton(unsigned r)
    : rank(r)
{ }

unsigned Skeleton::GetRank() const
{
    return rank;
}

Skeleton::~Skeleton() { }

//TODO: is correct num_cols and num_rows to mode=transpose????
/* calc is array wirh rank elem */
void Skeleton::MatvecF(double alpha, const double *vec,
        double *res, double *calc, mv_mode mode) const
{
    unsigned i, k, num_cols, num_rows;
    const Matrix& U = mode == normal ? GetU() : GetV();
    const Matrix& V = mode == normal ? GetV() : GetU();
    num_cols = GetNumCols();
    num_rows = GetNumRows();
    /* calc alpha * V.T * vec */
    for(k = 0; k < rank; ++k) {
        calc[k] = 0;
        for(i = 0; i < num_cols; ++i) {
            calc[k] += V.elem(i, k) * vec[i];
        }
        calc[k] *= alpha;
    }
    /* calc U * (alpha * V.T * vec) */
    for(i = 0; i < num_rows; ++i) {
        res[i] = 0;
        for(k = 0; k < rank; ++k) {
            res[i] += U.elem(i, k) * calc[k];
        }
    }
}

void Skeleton::MatvecUT(double alpha, const double *vec,
        double *res, mv_mode mode) const
{
    unsigned i, k, num_cols, num_rows;
    double tmp;
    const Matrix& U = mode == normal ? GetU() : GetV();
    const Matrix& V = mode == normal ? GetV() : GetU();
    num_cols = GetNumCols();
    num_rows = GetNumRows();
    for(i = 0; i < num_rows; ++i) {
        res[i] = 0.0;
    }
    for(k = 0; k < rank; ++k) {
        tmp = 0.0;
        for(i = 0; i < num_rows; ++i) {
            if(i < num_cols) {
                tmp += V.elem(num_rows - 1 - i, k) * vec[num_rows - 1 - i];
            }
            res[num_rows - 1 - i] += tmp * U.elem(num_rows - 1 - i, k);
        }
    }
    for(i = 0; i < num_rows; ++i) {
        res[i] *= alpha;
    }
}

void Skeleton::MatvecLT(double alpha, const double *vec,
        double *res, mv_mode mode) const
{
    unsigned i, k, num_cols, num_rows;
    double tmp;
    const Matrix& U = mode == normal ? GetU() : GetV();
    const Matrix& V = mode == normal ? GetV() : GetU();
    num_cols = GetNumCols();
    num_rows = GetNumRows();
    for(i = 0; i < num_rows; ++i) {
        res[i] = 0.0;
    }
    for(k = 0; k < rank; ++k) {
        tmp = 0.0;
        for(i = 0; i < num_rows; ++i) {
            if(i < num_cols) {
                tmp += V.elem(i, k) * vec[i];
            }
            res[i] += tmp * U.elem(i, k);
        }
    }
    for(i = 0; i < num_rows; ++i) {
        res[i] *= alpha;
    }
}
