#ifndef SKELETON_HPP
#define SKELETON_HPP

#include "matrix.hpp"
#include "indexlist.hpp"

/*!
    \file skeleton.hpp
    \brief Rank/skeleton decomposition base class

    The file defines the base class of rank/skeleton decomposition base class
    with function to fast calculate matvec with matrix or its upper/lower
    triangular part if its stored in skeleton format.
*/

//! Abstract class of skeleton/rank decomposition of matrix
/*! Class stored skeleton decomposition of real mxn matrix M
 * that is M = U * V.T where U and V is real matrices of sizes mxr and nxr.
 * This class implement functions to fast calculation of matvec
 * with M, M.T or their upper/lower triangular parts.
 */
class Skeleton : public Matrix {
protected:
    unsigned rank;
public:
    //! Initialize skeleton decompostion of rank r
    Skeleton(unsigned r);
    //! Get rank of matrix
    unsigned GetRank() const;
    //! Return matrix of which skeleton decomposition is made
    virtual const Matrix& GetM() const = 0;
    //! Return U factor of skeleton decomposition
    virtual const Matrix& GetU() const = 0;
    //! Return V factor of skeleton decomposition
    virtual const Matrix& GetV() const = 0;
    //! Matrix for which operations be used: matrix itself or its transpose
    enum mv_mode {normal, transpose};
    //! Matvec with Full part of matrix
    /*! Calculate res = alpha * A * vec,
     *  where A is mxn matrix (A equal M or M.T) with rank r.
     *  @param alpha real scalar multiplier
     *  @param vec double vector of size n
     *  @param res double vector of size m; the result will be save in res
     *  @param calc double vector of size r; required to calculation
     *  @param mode if mode=normal then A=M else A=M.T
    */
    void MatvecF(double alpha, const double *vec,
        double *res, double *calc, mv_mode mode = normal) const;
    //! Matvec with Upper Triangular part of matrix
    /*! Calculate res = alpha * upper_triangular_part(A) * vec,
     *  where A is mxn matrix (A equal M or M.T),
     *  upper_triangular_part(A) is a mxn matrix which upper triangular part
     *  equal to A.
     *  @param alpha real scalar multiplier
     *  @param vec double vector of size n
     *  @param res double vector of size m; the result will be save in res
     *  @param mode if mode=normal then A=M else A=M.T
    */
    void MatvecUT(double alpha, const double *vec,
        double *res, mv_mode mode = normal) const;
    //! Matvec with Lower Triangular part of matrix
    /*! Calculate res = alpha * lower_triangular_part(A) * vec,
     *  where A is mxn matrix (A equal M or M.T),
     *  lower_triangular_part(A) is a mxn matrix which lower triangular part
     *  equal to A.
     *  @param alpha real scalar multiplier
     *  @param vec double vector of size n
     *  @param res double vector of size m; the result will be save in res
     *  @param mode if mode=normal then A=M else A=M.T
    */
    void MatvecLT(double alpha, const double *vec,
        double *res, mv_mode mode = normal) const;
    //! The destructor
    ~Skeleton();
};

#endif /* SKELETON_HPP */
