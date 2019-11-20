#pragma once
#include <iostream>
#include <cmath>
#include "TMatrix.h"
#ifdef __CMAKE__
	#define _ASSERT(x) if(!(x)) {exit(-1);}
#endif // !_ASSERT

const long DOUBLE_IN_BYTES = sizeof(double) / sizeof(char);
const double pi = 3.1415926;


#define FreeIf(buffer) if(buffer){delete[] buffer; buffer = nullptr;}

template<typename T>
bool AbsGreater(T a, T b);

class TMatrix;
class TCircularBuffer {
private:
	double* data;
	int  m_nLength;
	int m_nCurrentIndex;
	void allocSpace(double dInitialValue = 0.0f);
public:
	TCircularBuffer(int nLength);
	TCircularBuffer(double* dat, int nLength);
	TCircularBuffer(const TCircularBuffer& buffer);
	TCircularBuffer(const TMatrix& mat);
	TCircularBuffer & push(double num);
	double last();
	int length() const;
	TCircularBuffer& operator+=(double value);
	TCircularBuffer& operator-=(double value);
	TCircularBuffer& operator*=(double value);
	TCircularBuffer operator * (double value);
	TCircularBuffer& operator/=(int value);
	TCircularBuffer& reverse();
	TCircularBuffer toReverse();
	TCircularBuffer& operator+=(const TCircularBuffer& buffer);
	TCircularBuffer operator-(); 
	TCircularBuffer slice(int start, int end) const;
	const double& operator[](int index) const;
	TMatrix toMatrix();
	TMatrix differ(bool bBackward = true);
	friend std::ostream& operator << (std::ostream& os, const TCircularBuffer& buffer);
	friend std::istream& operator>>(std::istream& in, TCircularBuffer& buffer);
	class const_iterator {
		int m_nIndex;
		TCircularBuffer& m_buffer;
	public:
		const_iterator(TCircularBuffer& buffer, int nIndex = 0) : m_buffer(buffer), m_nIndex(nIndex) { }
		bool operator<(const const_iterator& iter) {
			return m_nIndex < iter.m_nIndex;
		}
		bool operator<=(const const_iterator& iter) {
			return m_nIndex <= iter.m_nIndex;
		}
		const const_iterator& operator=(const const_iterator& iter) {
			m_nIndex = iter.m_nIndex;
			m_buffer = iter.m_buffer;
			return *this;
		}
		const_iterator& operator++() {
			++m_nIndex;
			return *this;
		}
		const_iterator& operator++(int) {
			++m_nIndex;
			return *this;
		}
		const double& operator*() {
			auto index = (m_buffer.m_nCurrentIndex + m_nIndex - 1 + m_buffer.m_nLength) % m_buffer.m_nLength;
			return m_buffer.data[index];
		}
	};
	const_iterator begin();
	const_iterator end();
};
TCircularBuffer operator+(const TCircularBuffer& bufferA, const TCircularBuffer& bufferB);

template<typename T>
inline bool AbsGreater(T a, T b)
{
	return abs(a) > abs(b);
}