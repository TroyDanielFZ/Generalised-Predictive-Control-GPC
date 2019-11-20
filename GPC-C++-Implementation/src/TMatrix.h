#ifndef __MATRIX_H__
#define __MATRIX_H__

// ----------------------------------------------------------------------------------
// Note before coding
// There are quite many ways to implement a matrix class, for example, 1D/2D arrays
// At the very beginning, I want to implement it with 1D array, since it's easy to 
// handle the memories. However, notice that, the basic operations of matrix are 
// inverse, which shall take loads of row-switch operations. In this sense, a 2D
// array seems work better, despite the cost for memory management.
// ----------------------------------------------------------------------------------
#include <iostream>
#include <ostream>

class TCircularBuffer;
class TMatrix {
private:
	int m_nRow, m_nCol;
	int m_nTotalSize;
	double** m_pData;

	void allocSpace(double dInitialValue = 0.0f);
	void reallocSpace(int nRow, int nCol);
	TMatrix expHelper(const TMatrix&, int);
public:
	friend class TCircularBuffer;
	// constructors, copy constructors and destructors
	TMatrix(int rows, int cols);
	TMatrix();
	TMatrix(double* dat, int rows, int cols);
	TMatrix(const TMatrix&);
	~TMatrix();
	TMatrix& operator=(const TMatrix&);

	// return the total rows
	const int cols() const;
	// return the total columns
	const int rows() const;

	// ------------------------------------
	// Memeber Asscessors
	// ------------------------------------
	// return the element reference at row x, column y
	inline double& operator()(int x, int y) { return m_pData[x][y]; }
	// return the pointer to the row r
	double* operator[](int r);

	// math operators
	TMatrix& operator+=(const TMatrix&);
	TMatrix& operator+=(double value);
	TMatrix& operator-=(const TMatrix&);
	TMatrix& operator-=(double value);
	TMatrix& operator*=(const TMatrix&);
	TMatrix& operator*=(double);
	TMatrix& operator/=(double);
	friend TMatrix operator-(double, const TMatrix&);
	TMatrix operator-();// unary minus operator
	TMatrix  operator^(int);



	friend std::ostream& operator<<(std::ostream&, const TMatrix&);
	friend std::istream& operator>>(std::istream&, TMatrix&);

	void swapRows(int, int);
	TMatrix transpose();

	static TMatrix createIdentity(int);

	// functions on augmented matrices
	// Maybe, it's better to make it private
	static TMatrix augment(TMatrix, TMatrix, bool bHorizontal = true);
	// Act gaussian elimination on this matrix
	TMatrix& gaussianEliminate();

	// return sub matrix of this matrix, support positive indice and/or 
	// negative indice, for example
	//		A = B.slice(1, 2, 0, 1)
	// will give the submatrix of row 1 and 2, column 0 and 1. Noting 
	// that, the rows and columns start from zero
	//		A = B.slice(-2, -1, -2, -1)
	// will give the submatrix from the last but one row to last row, 
	// last but one and last column. Noting that the negative indice 
	// is starting from -1.
	// A mixture of positive/negative indice is supported.
	TMatrix slice(int row_start, int row_end, int col_start, int col_end);
	// give the inverse of this matrix
	TMatrix inverse();
	// Given the solution $X$ for $A X = B$
	TMatrix divideBy(const TMatrix& denominator);

	// ------------------------------------
	// RowElement, for better row assignment
	// ------------------------------------
	class RowElement {
	private:
		int _index;
		TMatrix& _parent;
	public:
		// constructor
		RowElement(TMatrix& parent, int index);
		RowElement(const RowElement& rowElement);

		const int length() const{ return _parent.cols(); }
		double& operator[](int index);
		// operator= overide
		const RowElement& operator=(const RowElement& rowElement);
		const RowElement& operator=(const TCircularBuffer& buffer); 
		// convertor
		TCircularBuffer toCircularBuffer() const; 
		// IO operation
		friend std::ostream& operator<<(std::ostream& os, const RowElement& row);
	};
	RowElement Row(int row_index);
};

// binary math operators
TMatrix operator+(const TMatrix&, const TMatrix&);
TMatrix operator-(const TMatrix&, const TMatrix&);
TMatrix operator*(const TMatrix&, const TMatrix&);
TMatrix operator*(const TMatrix&, double);
TMatrix operator*(double, const TMatrix&);
TMatrix operator/(const TMatrix&, double);

#endif
