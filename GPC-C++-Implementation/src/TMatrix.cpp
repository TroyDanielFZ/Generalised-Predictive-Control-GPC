#include <stdexcept>
#include <iomanip>
#include "TMatrix.h"
#include "TMath.h"
#define EPS 1e-10

using std::ostream; 
using std::istream; 
using std::endl;
using std::domain_error;

TMatrix::TMatrix(int rows, int cols) 
	: m_nRow(rows)
	, m_nCol(cols)
	, m_nTotalSize(rows * cols)
	, m_pData(nullptr){
    allocSpace(0.0f);
}

TMatrix::TMatrix() : m_nRow(0),
	m_nCol(0),
	m_nTotalSize(0),
	m_pData(nullptr){
}

TMatrix::TMatrix(double* dat, int rows, int cols)
	: m_nRow(rows)
	, m_nCol(cols)
	, m_nTotalSize(rows * cols)
	, m_pData(nullptr)
{
    allocSpace(0.0f);
	int nRowBytes = sizeof(double) / sizeof(char) * m_nCol;
    for (int i = 0; i < m_nRow; ++i) {
		memcpy(m_pData[i], dat+i*cols, nRowBytes);
    }

}

TMatrix::~TMatrix() {
	if (m_pData) for (int i = 0; i < m_nRow; ++i) FreeIf(m_pData[i]);
	FreeIf(m_pData);
	m_nRow= m_nCol= m_nTotalSize = 0;
}

TMatrix::TMatrix(const TMatrix& m) : m_nRow(m.m_nRow), m_nCol(m.m_nCol) {
    allocSpace();
	int nRowBytes = sizeof(double) / sizeof(char) * m_nCol;
    for (int i = 0; i < m_nRow; ++i) {
		memcpy(m_pData[i], m.m_pData[i], nRowBytes);
    }
}

TMatrix& TMatrix::operator=(const TMatrix& m) {
    if (this == &m) {
        return *this;
    }

    if (m_nRow != m.m_nRow || m_nCol != m.m_nCol) {
		reallocSpace(m.m_nRow, m.m_nCol);
    }

	int nRowBytes = sizeof(double) / sizeof(char) * m_nCol;
    for (int i = 0; i < m_nRow; ++i) {
		memcpy(m_pData[i], m.m_pData[i], nRowBytes);
    }
    return *this;
}

const int TMatrix::cols() const
{
	return m_nCol;
}
const int TMatrix::rows() const
{
	return m_nRow;
}
TMatrix& TMatrix::operator+=(const TMatrix& m) {
    for (int i = 0; i < m_nRow; ++i) {
        for (int j = 0; j < m_nCol; ++j) {
            m_pData[i][j] += m.m_pData[i][j];
        }
    }
    return *this;
}

TMatrix& TMatrix::operator+=(double value) {
	for (int r = 0; r < m_nRow; ++r) {
		double* dat = m_pData[r];
		for (int c = 0; c < m_nCol; ++c) {
			dat[c] += value;
		}
	}
	return *this;
}

TMatrix& TMatrix::operator-=(const TMatrix& m) {
    for (int i = 0; i < m_nRow; ++i) {
        for (int j = 0; j < m_nCol; ++j) {
            m_pData[i][j] -= m.m_pData[i][j];
        }
    }
    return *this;
}

TMatrix& TMatrix::operator-=(double value) {
    for (int i = 0; i < m_nRow; ++i) {
        for (int j = 0; j < m_nCol; ++j) {
            m_pData[i][j] -= value;
        }
    }
    return *this;
	// TODO: insert return statement here
}

// TODO: reference count
TMatrix& TMatrix::operator*=(const TMatrix& m) {
	_ASSERT(m_nCol == m.m_nRow);
    TMatrix temp(m_nRow, m.m_nCol);
    for (int i = 0; i < temp.m_nRow; ++i) {
        for (int j = 0; j < temp.m_nCol; ++j) {
            for (int k = 0; k < m_nCol; ++k) {
                temp.m_pData[i][j] += (m_pData[i][k] * m.m_pData[k][j]);
            }
        }
    }
    return (*this = temp);
}

TMatrix& TMatrix::operator*=(double num)
{
    for (int i = 0; i < m_nRow; ++i) {
        for (int j = 0; j < m_nCol; ++j) {
            m_pData[i][j] *= num;
        }
    }
    return *this;
}

TMatrix& TMatrix::operator/=(double num) {
    for (int i = 0; i < m_nRow; ++i) {
        for (int j = 0; j < m_nCol; ++j) {
            m_pData[i][j] /= num;
        }
    }
    return *this;
}

TMatrix TMatrix::operator-() {
    TMatrix temp(m_nRow, m_nCol);
	for (int r = 0; r < m_nRow; ++r) {
		double* dat = temp[r];
		double* dat_src = m_pData[r];
		for (int c = 0; c < m_nCol; ++c) {
			dat[c] = -dat_src[c];
		}
	}
	return  temp;
}

TMatrix TMatrix::operator^(int num) {
    TMatrix temp(*this);
    return expHelper(temp, num);
}

double* TMatrix::operator[](int row) {
	return row < m_nRow ? m_pData[row] : nullptr;
}

void TMatrix::swapRows(int r1, int r2)
{
    double *temp = m_pData[r1];
    m_pData[r1] = m_pData[r2];
    m_pData[r2] = temp;
}

TMatrix TMatrix::transpose() {
    TMatrix ret(m_nCol, m_nRow);
    for (int i = 0; i < m_nRow; ++i) {
        for (int j = 0; j < m_nCol; ++j) {
            ret.m_pData[j][i] = m_pData[i][j];
        }
    }
    return ret;
}


/* STATIC CLASS FUNCTIONS
 ********************************/

TMatrix TMatrix::createIdentity(int size) {
    TMatrix temp(size, size);
	for (int i = 0; i < temp.m_nRow; ++i)
		temp.m_pData[i][i] = 1;
    return temp;
}

// TODO:
// it seems that there may be some improvements
//Matrix Matrix::solve(Matrix A, Matrix b) {
//    // Gaussian elimination
//    for (int i = 0; i < A.m_nRow; ++i) {
//        if (A.m_pData[i][i] == 0) { // pivot 0 - throw error
//            throw domain_error("Error: the coefficient matrix has 0 as a pivot. Please fix the input and try again.");
//        }
//        for (int j = i + 1; j < A.m_nRow; ++j) {
//            for (int k = i + 1; k < A.m_nCol; ++k) {
//                A.m_pData[j][k] -= A.m_pData[i][k] * (A.m_pData[j][i] / A.m_pData[i][i]);
//                if (A[j][k] < EPS && A[j][k] > -1*EPS)
//                    A[j][k] = 0;
//            }
//            b.m_pData[j][0] -= b.m_pData[i][0] * (A.m_pData[j][i] / A.m_pData[i][i]);
//            if (A[j][0] < EPS && A.m_pData[j][0] > -1*EPS)
//                A[j][0] = 0;
//            A[j][i] = 0;
//        }
//    }
//
//    // Back substitution
//    Matrix x(b.m_nRow, 1);
//    x.m_pData[x.m_nRow - 1][0] = b.m_pData[x.m_nRow - 1][0] / A.m_pData[x.m_nRow - 1][x.m_nRow - 1];
//    if (x.m_pData[x.m_nRow - 1][0] < EPS && x.m_pData[x.m_nRow - 1][0] > -1*EPS)
//        x.m_pData[x.m_nRow - 1][0] = 0;
//    for (int i = x.m_nRow - 2; i >= 0; --i) {
//        int sum = 0;
//        for (int j = i + 1; j < x.m_nRow; ++j) {
//            sum += A.m_pData[i][j] * x.m_pData[j][0];
//        }
//        x.m_pData[i][0] = (b.m_pData[i][0] - sum) / A.m_pData[i][i];
//        if (x.m_pData[i][0] < EPS && x.m_pData[i][0] > -1*EPS)
//            x.m_pData[i][0] = 0;
//    }
//
//    return x;
//}
//
//Matrix Matrix::bandSolve(Matrix A, Matrix b, int k) {
//    // optimized Gaussian elimination
//    int bandsBelow = (k - 1) / 2;
//    for (int i = 0; i < A.m_nRow; ++i) {
//        if (A.m_pData[i][i] == 0) { // pivot 0 - throw exception
//            throw domain_error("Error: the coefficient matrix has 0 as a pivot. Please fix the input and try again.");
//        }
//        for (int j = i + 1; j < A.m_nRow && j <= i + bandsBelow; ++j) {
//            int k = i + 1;
//            while (k < A.m_nCol && A.m_pData[j][k]) {
//                A.m_pData[j][k] -= A.m_pData[i][k] * (A.m_pData[j][i] / A.m_pData[i][i]);
//                k++;
//            }
//            b.m_pData[j][0] -= b.m_pData[i][0] * (A.m_pData[j][i] / A.m_pData[i][i]);
//            A.m_pData[j][i] = 0;
//        }
//    }
//
//    // Back substitution
//    Matrix x(b.m_nRow, 1);
//    x.m_pData[x.m_nRow - 1][0] = b.m_pData[x.m_nRow - 1][0] / A.m_pData[x.m_nRow - 1][x.m_nRow - 1];
//    for (int i = x.m_nRow - 2; i >= 0; --i) {
//        int sum = 0;
//        for (int j = i + 1; j < x.m_nRow; ++j) {
//            sum += A.m_pData[i][j] * x.m_pData[j][0];
//        }
//        x.m_pData[i][0] = (b.m_pData[i][0] - sum) / A.m_pData[i][i];
//    }
//
//    return x;
//}
//
// functions on VECTORS
//double Matrix::dotProduct(Matrix a, Matrix b) {
//    double sum = 0;
//    for (int i = 0; i < a.m_nRow; ++i) {
//        sum += (a(i, 0) * b(i, 0));
//    }
//    return sum;
//}
//
// functions on AUGMENTED matrices
TMatrix TMatrix::augment(TMatrix A, TMatrix B, bool bHorizontal) {
	if(bHorizontal){
		if (A.m_nRow != B.m_nRow) throw("The row of augmented matrice A and B must be the same.");
		TMatrix AB(A.m_nRow, A.m_nCol + B.m_nCol);
		for (int i = 0; i < AB.m_nRow; ++i) {
			memcpy(AB[i], A[i], DOUBLE_IN_BYTES * A.m_nCol);
			memcpy(AB[i] + A.m_nCol, B[i], DOUBLE_IN_BYTES * B.m_nCol);
		}
		return AB;
	}
	else {
		if (A.m_nCol != B.m_nCol) throw("The row of augmented matrice A and B must be the same.");
		TMatrix AB(A.m_nRow + B.m_nRow,  A.m_nCol);
		for (int i = 0; i < A.m_nRow; ++i) {
			memcpy(AB[i], A[i], DOUBLE_IN_BYTES * A.m_nCol);
		}
		for (int i = A.m_nRow; i < AB.m_nRow; ++i) {
			memcpy(AB[i], B[i - A.m_nRow], DOUBLE_IN_BYTES * B.m_nCol);
		}
		return AB;
	}
}

TMatrix & TMatrix::gaussianEliminate() {
	int rows = m_nRow, cols = m_nCol;
	if (rows > cols)
		throw("Gaussian Eliminate function is designed to "
			"perform Gaussian Elimination for inverse calculation."
			"The equation is $Ax=B$, and the elimnation was performed "
			"on $[A, B]$.");
	auto data = m_pData;
    int Acols = cols - 1;

    int row = 0; // row tracker

    // iterate through the rows to get the row-reduce form
	for (int row = 0; row < rows; ++row) {
        // find a pivot for the row
		int max_row = row;	// row index for pivot row
		double max_value = data[row][row];
		for (int i = row + 1; i < rows; ++i) { // seek for pivot row
			if (AbsGreater<double>(data[i][row], max_value)) {
				max_row = i;
				max_value = data[i][row];
			}
		}
		if (!AbsGreater<double>(max_value, EPS)) 
			throw("Gaussian elimination failed for non-full-rank matrix.");
		// switch the pivot to with current row
		if (max_row != row) swapRows(max_row, row);
        bool pivot_found = false;

		// divide the pivot row by the pivot value, e.g. max_value
		double* pivot_row = data[row];
		for (int col = row; col < cols; ++col) 
			pivot_row[col] = pivot_row[col] / max_value;
		// sub the rows below the pivot row
		for (int i = row + 1; i < rows; ++i) {
			double* cur_row = data[i];
			if (!AbsGreater<double>(cur_row[row], EPS)) continue; // beyond the precision
			// perform elimination as normal if pivot was found
			double ratio = cur_row[row];
			for (int col = row; col < cols; ++col)
				cur_row[col] -= pivot_row[col] * ratio;
		}
	}
	// change the row-reduce form into [I Inverse] form
	for (int row = rows - 1; row > 0; --row) {
		double* pivot_row = data[row];
		for (int i = row - 1; i >= 0; --i) {
			double* cur_row = data[i];
			if (!AbsGreater<double>(cur_row[row], EPS)) continue; // beyond the precision
			// perform elimination as normal if pivot was found
			for(int col = cols - 1; col >= row; -- col)
			//for (int col = row; col < cols; ++col)
				cur_row[col] -=  pivot_row[col] * cur_row[row];
		}

	}
	return *this;
}
// If the index is negative, it countdown from the ends
TMatrix TMatrix::slice(int row_start, int row_end, int col_start, int col_end) {
	if (row_start < 0) row_start += m_nRow;
	if (row_end < 0) row_end += m_nRow;
	if (col_start < 0) col_start += m_nCol;
	if (col_end < 0) col_end += m_nCol;
	TMatrix temp(row_end - row_start + 1, col_end - col_start + 1);
	int nDoubleInBytes = sizeof(double) / sizeof(char);
	for (int i = row_start; i <= row_end; ++i) {
		memcpy(temp[i - row_start], (*this)[i] + col_start, nDoubleInBytes * (col_end - col_start + 1));
	}
	return temp;
}
TMatrix TMatrix::inverse() {
	if (m_nRow != m_nRow || m_nRow <= 0) throw("Non-square matrix or empty matrix doesn't have an inverse.");
	return augment(*this, createIdentity(m_nRow)).gaussianEliminate().slice(0, m_nRow-1, m_nCol, 2 * m_nCol -1);
}
TMatrix TMatrix::divideBy(const TMatrix& B) {
	if (m_nRow != m_nRow || m_nRow <= 0) throw("Non-square matrix or empty matrix isn't inversible.");
	//return augment(*this, B).gaussianEliminate().slice(0, m_nRow-1, m_nCol, m_nCol + B.m_nCol -1);
	return augment(*this, B).gaussianEliminate().slice(0, -1, m_nCol, -1);
}
TMatrix::RowElement TMatrix::Row(int row_index)
{
	return RowElement(*this, row_index);
}
//
//Matrix Matrix::rowReduceFromGaussian() {
//    Matrix R(*this);
//    int rows = R.m_nRow;
//    int cols = R.m_nCol;
//
//    int i = rows - 1; // row tracker
//    int j = cols - 2; // column tracker
//
//    // iterate through every row
//    while (i >= 0)
//    {
//        // find the pivot column
//        int k = j - 1;
//        while (k >= 0) {
//            if (R(i, k) != 0)
//                j = k;
//            k--;
//        }
//
//        // zero out elements above pivots if pivot not 0
//        if (R(i, j) != 0) {
//       
//            for (int t = i - 1; t >= 0; --t) {
//                for (int s = 0; s < cols; ++s) {
//                    if (s != j) {
//                        R(t, s) = R(t, s) - R(i, s) * (R(t, j) / R(i, j));
//                        if (R(t, s) < EPS && R(t, s) > -1*EPS)
//                            R(t, s) = 0;
//                    }
//                }
//                R(t, j) = 0;
//            }
//
//            // divide row by pivot
//            for (int k = j + 1; k < cols; ++k) {
//                R(i, k) = R(i, k) / R(i, j);
//                if (R(i, k) < EPS && R(i, k) > -1*EPS)
//                    R(i, k) = 0;
//            }
//            R(i, j) = 1;
//
//        }
//
//        i--;
//        j--;
//    }
//
//    return R;
//}
//
//void Matrix::readSolutionsFromRREF(ostream& os) {
//    Matrix R(*this);
//
//    // print number of solutions
//    bool hasSolutions = true;
//    bool doneSearching = false;
//    int i = 0;
//    while (!doneSearching && i < m_nRow) {
//        bool allZeros = true;
//        for (int j = 0; j < m_nCol - 1; ++j) {
//            if (R(i, j) != 0)
//                allZeros = false;
//        }
//        if (allZeros && R(i, m_nCol - 1) != 0) {
//            hasSolutions = false;
//            os << "NO SOLUTIONS" << endl << endl;
//            doneSearching = true;
//        } else if (allZeros && R(i, m_nCol - 1) == 0) {
//            os << "INFINITE SOLUTIONS" << endl << endl;
//            doneSearching = true;
//        } else if (m_nRow < m_nCol - 1) {
//            os << "INFINITE SOLUTIONS" << endl << endl;
//            doneSearching = true;
//        }
//        i++;
//    }
//    if (!doneSearching)
//        os << "UNIQUE SOLUTION" << endl << endl;
//
//    // get solutions if they exist
//    if (hasSolutions) {
//        Matrix particular(m_nCol - 1, 1);
//        Matrix special(m_nCol - 1, 1);
//
//        for (int i = 0; i < m_nRow; ++i) {
//            bool pivotFound = false;
//            bool specialCreated = false;
//            for (int j = 0; j < m_nCol - 1; ++j) {
//                if (R(i, j) != 0) {
//                    // if pivot variable, add b to particular
//                    if (!pivotFound) {
//                        pivotFound = true;
//                        particular(j, 0) = R(i, m_nCol - 1);
//                    } else { // otherwise, add to special solution
//                        if (!specialCreated) {
//                            special = Matrix(m_nCol - 1, 1);
//                            specialCreated = true;
//                        }
//                        special(j, 0) = -1 * R(i, j);
//                    }
//                }
//            }
//            os << "Special solution:" << endl << special << endl;
//        }
//        os << "Particular solution:" << endl << particular << endl;
//    }
//}
//
//Matrix Matrix::inverse() {
//    Matrix I = Matrix::createIdentity(m_nRow);
//    Matrix AI = Matrix::augment(*this, I);
//    Matrix U = AI.gaussianEliminate();
//    Matrix IAInverse = U.rowReduceFromGaussian();
//    Matrix AInverse(m_nRow, m_nCol);
//    for (int i = 0; i < AInverse.m_nRow; ++i) {
//        for (int j = 0; j < AInverse.m_nCol; ++j) {
//            AInverse(i, j) = IAInverse(i, j + m_nCol);
//        }
//    }
//    return AInverse;
//}
//

/* PRIVATE HELPER FUNCTIONS
 ********************************/

void TMatrix::allocSpace(double dInitialValue) {
    m_pData = new double*[m_nRow];
    for (int i = 0; i < m_nRow; ++i) {
		m_pData[i] = new double[m_nCol] {dInitialValue};
    }
}

void TMatrix::reallocSpace(int nRow, int nCol) {
	if (m_pData) for (int i = 0; i < m_nRow; ++i) FreeIf(m_pData[i]);
	FreeIf(m_pData);
	m_nRow = nRow;
	m_nCol = nCol;
	m_nTotalSize = nRow * nCol;
	allocSpace();
}

TMatrix TMatrix::expHelper(const TMatrix& m, int num) {
    if (num == 0) { 
        return createIdentity(m.m_nRow);
    } else if (num == 1) {
        return m;
    } else if (num % 2 == 0) {  // num is even
        return expHelper(m * m, num/2);
    } else {                    // num is odd
        return m * expHelper(m * m, (num-1)/2);
    }
}

/* NON-MEMBER FUNCTIONS
 ********************************/

std::ostream& operator<<(std::ostream& os, const TMatrix::RowElement& row) {
	double* dat = row._parent[row._index];
	for (int i = 0; i < row.length(); ++i)
		std::cout << dat[i] << " ";
	std::cout << std::endl;
	return os;
}

TMatrix operator+(const TMatrix& m1, const TMatrix& m2) {
    TMatrix temp(m1);
    return (temp += m2);
}

TMatrix operator-(const TMatrix& m1, const TMatrix& m2) {
    TMatrix temp(m1);
    return (temp -= m2);
}

TMatrix operator-(double value, const TMatrix& m) {
    TMatrix temp(m);
	for(int r = 0; r < m.m_nRow; ++r){
		double* dat = temp[r];
		for (int c = 0; c < m.m_nCol; ++c)
			dat[c] = value - dat[c];
	}
	return temp;
}

TMatrix operator*(const TMatrix& m1, const TMatrix& m2) {
    TMatrix temp(m1);
    return (temp *= m2);
}

TMatrix operator*(const TMatrix& m, double num) {
    TMatrix temp(m);
    return (temp *= num);
}

TMatrix operator*(double num, const TMatrix& m) {
    return (m * num);
}

TMatrix operator/(const TMatrix& m, double num) {
    TMatrix temp(m);
    return (temp /= num);
}

ostream& operator<<(ostream& os, const TMatrix& m) {
    for (int i = 0; i < m.m_nRow; ++i) {
        os << std::right << std::setw(12) << std::setprecision(4)<< m.m_pData[i][0];
        for (int j = 1; j < m.m_nCol; ++j) {
            os << " " << std::right << std::setw(8) << std::setprecision(4) << m.m_pData[i][j];
        }
        os << endl;
    }
    return os;
}

istream& operator>>(istream& is, TMatrix& m) {
    for (int i = 0; i < m.m_nRow; ++i) {
        for (int j = 0; j < m.m_nCol; ++j) {
            is >> m.m_pData[i][j];
        }
    }
    return is;
}

TMatrix::RowElement::RowElement(TMatrix& parent, int index)
	: _parent(parent), _index(index) { }
TMatrix::RowElement::RowElement(const RowElement& rowElement)
	:_index(rowElement._index), _parent(rowElement._parent) { }
const TMatrix::RowElement& TMatrix::RowElement::operator=(const RowElement& rowElement)
{
	if (_parent.m_nCol != rowElement.length()) throw "The size doesn't match!";
	_parent = rowElement._parent;
	double* dat_src = rowElement._parent.m_pData[rowElement._index];
	double* dat = _parent.m_pData[_index];
	for (int i = 0; i < _parent.m_nCol; ++i)
		dat[i] = dat_src[i];
	return *this;
}

const TMatrix::RowElement& TMatrix::RowElement::operator=(const TCircularBuffer& buffer) {
	if (_parent.m_nCol != buffer.length()) throw "The size doesn't match!";
	double* dat = _parent.m_pData[_index];
	for (int i = 0; i < _parent.m_nCol; ++i) {
		dat[i] = buffer[i];
	}
	return *this;
}

double& TMatrix::RowElement::operator[](int index)
{
	return _parent.m_pData[_index][index];
}

TCircularBuffer TMatrix::RowElement::toCircularBuffer() const {
	TCircularBuffer buffer(_parent.m_pData[_index], length());
	return buffer;
}
