#ifdef __THIS_FILE_ARE_NOT_USED__
#include "TVector.h"
#include "TMath.h"

void TColumnVector::reallocSpace(int nLength, double dInitialValue) {
	FreeIf(m_pData);
	m_pData = new double[nLength] {dInitialValue};
	m_nLength = nLength;
}

TColumnVector::TColumnVector(): m_pData(nullptr), m_nLength(0)
{
}

TColumnVector::TColumnVector(double* data, int nLenth): m_nLength(nLenth),m_pData(nullptr) {
	reallocSpace(nLenth);
	memcpy(m_pData, data, DOUBLE_IN_BYTES * nLenth);
}

double& TColumnVector::operator[](int index)
{
	return m_pData[index];
}

TColumnVector& TColumnVector::operator=(const TColumnVector& vecSrc) {
	if (m_nLength != vecSrc.m_nLength)  reallocSpace(vecSrc.m_nLength);
	memcpy(m_pData, vecSrc.m_pData, DOUBLE_IN_BYTES * m_nLength);
	return *this;
}

TColumnVector& TColumnVector::operator=(double value) {
	if (m_nLength < 1) reallocSpace(1, value);
	else for (int i = 0; i < m_nLength; ++i) m_pData[i] = value;
	return *this;
}

double TColumnVector::norm2() {
	double result = 0.0f;
	for (int i = 0; i < m_nLength; ++i)
		result += pow(m_pData[i], 2);
	return result;
}

double TColumnVector::norm_infinity() {
	if (m_nLength < 1) throw("Null vector doesn't have an infity norm");
	double result = abs(m_pData[0]);
	for (int i = 0; i < m_nLength; ++i) {
		double abs_value = abs(m_pData[i]);
		result = result > abs_value ? result : abs_value;
	}
	return result;
}

double TColumnVector::norm1() {
	double result = 0.0f;
	for (int i = 0; i < m_nLength; ++i)
		result += fabs(m_pData[i]);
	return result;
}

TColumnVector::~TColumnVector() {
	FreeIf(m_pData);
	m_nLength = 0;
}

std::ostream& operator<<(std::ostream& os, const TColumnVector& vec) {
	for (int i = 0; i < vec.m_nLength; ++i) {
		os << vec.m_pData[i] << std::endl;
	}
	return os;
}

std::istream& operator>>(std::istream& in, TColumnVector& vec) {
	for (int i = 0; i < vec.m_nLength; ++i)
		in >> vec.m_pData[i];
	return in;
}
#endif
