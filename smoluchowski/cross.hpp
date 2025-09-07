#ifndef CROSS_HPP
#define CROSS_HPP

#include "matrix.hpp"
#include "skeleton.hpp"
#include "indexlist.hpp"

/*!
    \file cross.hpp
    \brief Cross/CGR/CUR approximation

    The file defines Cross class which calculate skeleton
    decomposition of matrix by Cross approximation algorithm describe
    for example in
    "A Parallel Implementation of the Matrix Cross Approximation Method"
    D. A. Zheltkov and E. E. Tyrtyshnikov, 2015, https://www.mathnet.ru/vmp548.
*/



//TODO: add set max_rank to allocate U, V and marked_index with this size
//! Class of Cross approximation method to find skeleton decomposition of matrix
/*! Class find skeleton decomposition of real mxn matrix A
 * by adaptive Cross approximation method with intermidiate 1-rank
 * approximation.
 * Implementation based on the algorithm described in
 * "A Parallel Implementation of the Matrix Cross Approximation Method"
 * D. A. Zheltkov and E. E. Tyrtyshnikov, 2015, https://www.mathnet.ru/vmp548.
 */
class Cross : public Skeleton {
    const Matrix& A;
    ReallocatableMatrix U;
    ReallocatableMatrix V;
    double fnorm2; /*||U * V.T||_F^2*/
    double eps;

    /* indexes from approximation process */
    IndexList row_marked_index;
    IndexList col_marked_index;

    /* recalculate frobenius norm after add one
    * column to u and v */
    void RecalcFnorm();
    
public:
    //! Initialize Cross approximation for matrix M
    Cross(const Matrix& M);
    //! Find Cross approximation with tolerance paramater eps
    /* Calculate Cross approximation of mxn matrix M with stop criteria
     * eps * fnorm(U * V) >= |M_{ij}| * sqrt((m - rank) * (n - rank)).
     * \note Method may be called several times with different eps
     */
    void Approximate(double _eps);
    //! Get last used tolerance eps
    double GetEps() const;
    unsigned GetNumRows() const;
    unsigned GetNumCols() const;
    double elem(unsigned i, unsigned j) const;
    //! Get the original matrix
    const Matrix& GetM() const;
    const Matrix& GetU() const;
    const Matrix& GetV() const;
    //! The destructor
    ~Cross();
};

#endif /* CROSS_HPP */
