#ifndef MATRIX_HPP
#define MATRIX_HPP

#ifndef HANDLEXCEPTION
#define HANDLEXCEPTION 0
#endif

class Vector {
public:
    Vector() { }
    virtual ~Vector() { }
    virtual double elem(unsigned i) const = 0;
    virtual unsigned GetSize() const = 0;
};

class Matrix {
public:
    Matrix() { }
    virtual ~Matrix() { }
    virtual unsigned GetNumRows() const = 0;
    virtual unsigned GetNumCols() const = 0;
    virtual double elem(unsigned i, unsigned j) const = 0;

    class Row : public Vector {
        friend class Matrix;
        const Matrix* M;
        const unsigned i;
        Row(const Matrix* _M, unsigned _i);
    public:
        double elem(unsigned j) const;
        unsigned GetSize() const;
    };
    Row GetRow(unsigned i) const;
    class Column : public Vector {
        friend class Matrix;
        const Matrix* M;
        const unsigned j;
        Column(const Matrix* _M, unsigned _j);
    public:
        double elem(unsigned i) const;
        unsigned GetSize() const;
    };
    Column GetCol(unsigned j) const;
};


// binary expression with two matrix with equal size
class MatrixBinExpr : public Matrix {
protected:
    const Matrix& M1;
    const Matrix& M2;
public:
    MatrixBinExpr(const Matrix& _M1, const Matrix& _M2);
    // Matrix Calc() const; write calculated expr in matrix
    const Matrix& GetFirst() const;
    const Matrix& GetSecond() const;
    unsigned GetNumRows() const;
    unsigned GetNumCols() const;
};

class MatrixSum : public MatrixBinExpr {
public:
    MatrixSum(const Matrix& _M1, const Matrix& _M2);
    double elem(unsigned i, unsigned j) const;
};

class MatrixSub : public MatrixBinExpr {
public:
    MatrixSub(const Matrix& _M1, const Matrix& _M2);
    double elem(unsigned i, unsigned j) const;
};


enum {ReallocatableMatrix_add_allocsize = 10};

/* matrix stored by column so you can easy add columns*/
class ReallocatableMatrix : public Matrix {
    unsigned num_rows;
    unsigned num_cols;
    unsigned allocsize;
    double **mat;
public:
    ReallocatableMatrix(unsigned rows, unsigned start_cols);
    unsigned GetNumRows() const;
    unsigned GetNumCols() const;
    void Resize(unsigned needsize);
    double& elem(unsigned i, unsigned j);
    double elem(unsigned i, unsigned j) const;
    ~ReallocatableMatrix();
};

/* matrix with elems f(i,j)*/
class FunctionMatrix : public Matrix {
    unsigned num_rows;
    unsigned num_cols;
    double h;
    double (*f)(unsigned, unsigned, double);
public:
    FunctionMatrix(unsigned rows, unsigned cols, double _h,
            double (*func)(unsigned, unsigned, double));
    unsigned GetNumRows() const;
    unsigned GetNumCols() const;
    double elem(unsigned i, unsigned j) const;
};

#endif /* MATRIX_HPP */
