#ifndef SKELETON_HPP
#define SKELETON_HPP

#include "matrix.hpp"
#include "indexlist.hpp"

/* skeleton/rank decomposition of rank r:
 * U * V.T */
class Skeleton : public Matrix {
protected:
    unsigned rank;
public:
    Skeleton(unsigned r);
    unsigned GetRank() const;
    /* return original matrix or itself */
    virtual const Matrix& GetM() const = 0;
    virtual const Matrix& GetU() const = 0;
    virtual const Matrix& GetV() const = 0;
    enum mv_mode {normal, transpose};
    /* matvec with Full part of matrix U * V.T */
    void MatvecF(double alpha, const double *vec,
        double *res, double *calc, mv_mode mode = normal) const;
    /* matvec with Upper Triangular part */
    void MatvecUT(double alpha, const double *vec,
        double *res, mv_mode mode = normal) const;
    /* matvec with Lower Triangular part */
    void MatvecLT(double alpha, const double *vec,
        double *res, mv_mode mode = normal) const;
    ~Skeleton();
};

#endif /* SKELETON_HPP */
