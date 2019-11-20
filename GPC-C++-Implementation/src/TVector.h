#ifdef __THIS_FILE_ARE_NOT_USED__
#pragma once
#include <iostream>
class TMatrix;
typedef double TScalar;
class TRowVector;
class TColumnVector {
private:
	double* m_pData;
	int m_nLength;
	void reallocSpace(int nLength, double dInitialValue = 0.0f);
public:
	TColumnVector();
	TColumnVector(double* data, int nLenth);

	double & operator[](int index);
	TColumnVector& operator=(const TColumnVector& vecSrc);
	TColumnVector& operator=(double value);

	double norm2();
	double norm_infinity();
	double norm1();
	// maybe the following function will be implemented in the future
	//double norm_p(double p);

	friend std::ostream& operator<<(std::ostream&, const TColumnVector&);
	friend std::istream& operator>>(std::istream&, TColumnVector&);
	~TColumnVector();
};

#endif
