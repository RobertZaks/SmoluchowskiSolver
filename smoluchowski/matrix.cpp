#include "matrix.hpp"

/*=============Matrix===============*/

Matrix::Row::Row(const Matrix* _M, unsigned _i)
    : M(_M), i(_i)
{ }

double Matrix::Row::elem(unsigned j) const
{
    return M->elem(i, j);
}

unsigned Matrix::Row::GetSize() const
{
    return M->GetNumCols();
}

Matrix::Row Matrix::GetRow(unsigned i) const
{
    return Row(this, i);
}


Matrix::Column::Column(const Matrix* _M, unsigned _j)
    : M(_M), j(_j)
{ }

double Matrix::Column::elem(unsigned i) const
{
    return M->elem(i, j);
}

unsigned Matrix::Column::GetSize() const
{
    return M->GetNumRows();
}

Matrix::Column Matrix::GetCol(unsigned j) const
{
    return Column(this, j);
}


/*=============MatrixBinExpr===============*/

MatrixBinExpr::MatrixBinExpr(const Matrix& _M1, const Matrix& _M2)
    : M1(_M1), M2(_M2)
{ 
#if HANDLEXCEPTION
    if((M1.GetNumRows() != M2.GetNumRows()) || 
        (M1.GetNumCols() != M2.GetNumCols()))
    {
        throw "assume M1 and M2 have the same size";
    }
#endif
}

unsigned MatrixBinExpr::GetNumRows() const
{
    return M1.GetNumRows();
}

unsigned MatrixBinExpr::GetNumCols() const
{
    return M1.GetNumCols();
}

const Matrix& MatrixBinExpr::GetFirst() const
{
    return M1;
}

const Matrix& MatrixBinExpr::GetSecond() const
{
    return M2;
}

MatrixSum::MatrixSum(const Matrix& _M1, const Matrix& _M2)
    : MatrixBinExpr(_M1, _M2)
{ }

double MatrixSum::elem(unsigned i, unsigned j) const
{
#if HANDLEXCEPTION
    if((i >= GetNumRows()) || (j >= GetNumCols()))
        throw "assume i < num_rows and j < num_cols";
    }
#endif
    return M1.elem(i, j) + M2.elem(i, j);
}

MatrixSub::MatrixSub(const Matrix& _M1, const Matrix& _M2)
    : MatrixBinExpr(_M1, _M2)
{ }

double MatrixSub::elem(unsigned i, unsigned j) const
{
#if HANDLEXCEPTION
    if((i >= GetNumRows()) || (j >= GetNumCols()))
        throw "assume i < num_rows and j < num_cols";
    }
#endif
    return M1.elem(i, j) - M2.elem(i, j);
}



/*=============ReallocatableMatrix===============*/

ReallocatableMatrix::ReallocatableMatrix(unsigned rows, unsigned start_cols)
    : num_rows(rows), num_cols(start_cols), allocsize(start_cols), mat(0)
{
    Resize(allocsize);
}

void ReallocatableMatrix::Resize(unsigned needsize)
{
    unsigned i;
    double **tmp;
    if(needsize > allocsize) {
        allocsize = needsize + ReallocatableMatrix_add_allocsize;
        tmp = new double*[allocsize];
        for(i = 0; i < num_cols; ++i){
            tmp[i] = mat[i];
        }
        delete [] mat;
        mat = tmp;
    }
    for(i = num_cols; i < needsize; ++i) {
        mat[i] = new double[num_rows];
    }
    num_cols = needsize;
}

unsigned ReallocatableMatrix::GetNumRows() const
{
    return num_rows;
}

unsigned ReallocatableMatrix::GetNumCols() const
{
    return num_cols;
}

double& ReallocatableMatrix::elem(unsigned i, unsigned j)
{
#if HANDLEXCEPTION
    if((i < num_rows) && (j < num_cols))
        throw "assume i < num_rows and j < num_cols";
    }
#endif
    return mat[j][i];
}

double ReallocatableMatrix::elem(unsigned i, unsigned j) const
{
#if HANDLEXCEPTION
    if((i < num_rows) && (j < num_cols))
        throw "assume i < num_rows and j < num_cols";
    }
#endif
    return mat[j][i];
}

ReallocatableMatrix::~ReallocatableMatrix()
{
    unsigned i;
    for(i = 0; i < num_cols; ++i) {
        delete [] mat[i];
    }
    delete [] mat;
}


/*=============FunctionMatrix===============*/


FunctionMatrix::FunctionMatrix(unsigned rows, unsigned cols, double _h,
        double (*func)(unsigned, unsigned, double))
    : num_rows(rows), num_cols(cols), h(_h), f(func)
{ }


unsigned FunctionMatrix::GetNumRows() const
{
    return num_rows;
}

unsigned FunctionMatrix::GetNumCols() const
{
    return num_cols;
}

double FunctionMatrix::elem(unsigned i, unsigned j) const
{
#if HANDLEXCEPTION
    if((i < num_rows) && (j < num_cols))
        throw "assume i < num_rows and j < num_cols";
    }
#endif
    return f(i, j, h);
}
