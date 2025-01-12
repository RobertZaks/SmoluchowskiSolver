#include "cross.hpp"
#include <math.h>

Cross::Cross(const Matrix &M)
    : Skeleton(0), A(M), U(M.GetNumRows(), rank), 
    V(M.GetNumCols(), rank), fnorm2(0)
{ }

Cross::~Cross() { }

const Matrix& Cross::GetM() const
{
    return A;
}

const Matrix& Cross::GetU() const
{
    return U;
}

const Matrix& Cross::GetV() const
{
    return V;
}

unsigned Cross::GetNumRows() const
{
    return A.GetNumRows();
}

unsigned Cross::GetNumCols() const
{
    
    return A.GetNumCols();
}

double Cross::GetEps() const
{
    return eps;
}

double Cross::elem(unsigned i, unsigned j) const
{
    unsigned k;
    double res = 0;
#if HANDLEXCEPTION
    if((i < GetNumRows()) && (j < GetNumCols()))
        throw "assume i < num_rows and j < num_cols";
    }
#endif
    // calc (U * V.T)[i][j]
    for(k = 0; k < rank; ++k) {
        res += U.elem(i, k) * V.elem(j, k);
    }
    return res;
}


/* return argmax of abs value in Vector
 * ignoring some index */
static unsigned
argmax_ignore_idx(const Vector& vec, const IndexList& ignore_idx)
{
    unsigned res, i, j, size;
    double tmp, max_a = 0;
    size = ignore_idx.GetCurrentSize();
    i = ignore_idx.GetFirstNotListed();
    res = i;
    for(j = i; j < size; ++j) {
        for(; i < ignore_idx[j]; ++i) {
            tmp = fabs(vec.elem(i));
            if(tmp > max_a) {
                max_a = tmp;
                res = i;
            }
        }
        i = ignore_idx[j] + 1;
    }
    return res;
}

/* recalculate frobenius norm after add one
 * column to u and v */
void Cross::RecalcFnorm()
{
    unsigned i, k, num_cols, num_rows;
    double tmp1, tmp2;
    num_cols = GetNumCols();
    num_rows = GetNumRows();
    /* add dotproduts*/
    for(k = 0; k + 1 < rank; ++k) {
        for(tmp1 = 0, i = 0; i < num_rows; ++i) {
            tmp1 += U.elem(i, k) * U.elem(i, rank - 1);
        }
        for(tmp2 = 0, i = 0; i < num_cols; ++i) {
            tmp2 += V.elem(i, k) * V.elem(i, rank - 1);
        }
        fnorm2 += 2 * tmp1 * tmp2;
    }
    for(tmp1 = 0, i = 0; i < num_rows; ++i) {
        tmp1 += U.elem(i, rank - 1) * U.elem(i, rank - 1);
    }
    for(tmp2 = 0, i = 0; i < num_cols; ++i) {
        tmp2 += V.elem(i, rank - 1) * V.elem(i, rank - 1);
    }
    /* add ||u[r-1]||_2^2 * ||v[r-1]||_2^2 */
    fnorm2 += tmp1 * tmp2;
}

void Cross::Approximate(double _eps)
{
    unsigned i, j, t, num_cols, num_rows;
    double tmp;
    MatrixSub M(A, *this);
    num_cols = GetNumCols();
    num_rows = GetNumRows();
    eps = _eps;
    row_marked_index.AddIndex(num_rows);
    col_marked_index.AddIndex(num_cols);
    for(;;) {
        /* find argmax */
        i = argmax_ignore_idx(
                M.GetCol(col_marked_index.GetFirstNotListed()), 
                row_marked_index
            );
        j = argmax_ignore_idx(M.GetRow(i), col_marked_index);

        /* check criteria */
        tmp = M.elem(i, j);
        if(fabs(tmp) * 
            sqrt((num_rows - (rank + 1)) * (num_cols - (rank + 1))) <=
            eps * sqrt(fnorm2))
        {
            break;
        }

        /* add new columns to U and V */
        U.Resize(rank + 1);
        for(t = 0; t < num_rows; ++t) {
            U.elem(t, rank) = (M.elem(t, j)) / sqrt(fabs(tmp));
        }
        V.Resize(rank + 1);
        for(t = 0; t < num_cols; ++t) {
            V.elem(t, rank) = (M.elem(i, t)) * sqrt(fabs(tmp)) / tmp;
        }

        /* increase rank */
        rank++;

        /* recalc Frobenius norm*/
        RecalcFnorm();

        /* update marked_index */
        row_marked_index.AddIndex(i);
        col_marked_index.AddIndex(j);
    }
}

