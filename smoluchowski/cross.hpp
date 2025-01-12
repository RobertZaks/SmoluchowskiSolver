#ifndef CROSS_HPP
#define CROSS_HPP

#include "matrix.hpp"
#include "skeleton.hpp"
#include "indexlist.hpp"

/*TODO: add set max_rank to 
 * allocate U, V and marked_index with this size
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
    Cross(const Matrix& M);
    /* may be called several times with different _eps */
    void Approximate(double _eps);
    double GetEps() const;
    unsigned GetNumRows() const;
    unsigned GetNumCols() const;
    double elem(unsigned i, unsigned j) const;
    /* return original matrix */
    const Matrix& GetM() const;
    const Matrix& GetU() const;
    const Matrix& GetV() const;
    ~Cross();
};

#endif /* CROSS_HPP */
