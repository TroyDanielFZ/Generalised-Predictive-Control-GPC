#include "TMath.h"
std::ostream& operator<<(std::ostream& os, const TCircularBuffer& buffer)
{
	int splitter = (buffer.m_nCurrentIndex + 1) % buffer.m_nLength;
	os << buffer.m_nLength << ", " << buffer.m_nCurrentIndex << ": ";
	for (int i = buffer.m_nCurrentIndex; i < buffer.m_nLength; ++i)
		os << buffer.data[i] << ", ";
	for (int i = 0; i < buffer.m_nCurrentIndex; ++i)
		os << buffer.data[i] << ", ";
	return os;
}
std::istream& operator>>(std::istream& in, TCircularBuffer& buffer)
{
	char ch;
	in >> buffer.m_nLength >> ch >> buffer.m_nCurrentIndex >> ch;
	buffer.allocSpace();
	for (int i = buffer.m_nCurrentIndex; i < buffer.m_nLength; ++i)
		in >> buffer.data[i] >> ch;
	for (int i = 0; i < buffer.m_nCurrentIndex; ++i)
		in >> buffer.data[i] >> ch;
	return in;
}
TCircularBuffer operator+(const TCircularBuffer& bufferA, const TCircularBuffer& bufferB) {
	TCircularBuffer temp(bufferA);
	return temp += bufferB;
}

void TCircularBuffer::allocSpace(double dInitialValue)
{
	FreeIf(data);
	m_nCurrentIndex = 0;
	data = new double[m_nLength] {dInitialValue};
}

TCircularBuffer::TCircularBuffer(int nLength)
	:data(nullptr), m_nLength(nLength) {
	allocSpace();
}

TCircularBuffer::TCircularBuffer(double* dat, int nLength)
	:data(nullptr), m_nLength(nLength) {
	allocSpace();
	for (int i = 0; i < nLength; ++i) data[i] = dat[i];
}
TCircularBuffer::TCircularBuffer(const TCircularBuffer& buffer)
	:m_nLength(buffer.m_nLength)
	, data(new double[buffer.m_nLength])
	, m_nCurrentIndex(buffer.m_nCurrentIndex) {
	for (int i = 0; i < m_nLength; ++i)
		data[i] = buffer.data[i];
}

TCircularBuffer::TCircularBuffer(const TMatrix& mat):data(nullptr)
{
	TMatrix& m = const_cast<TMatrix&>(mat); // I know, the mat is const, and I won't change the value
	m_nLength = mat.m_nCol == 1 ? mat.m_nRow : mat.m_nCol;
	m_nCurrentIndex = 0;
	allocSpace();
	if (mat.m_nCol == 1) {
		for (int i = 0; i < m_nLength; ++i) {
			data[i] = m[i][0];
		}
	}
	else {
		double* dat = m[0];
		for (int i = 0; i < m_nLength; ++i) {
			data[i] = dat[i];
		}
	}
}

TCircularBuffer& TCircularBuffer::push(double num)
{
	data[m_nCurrentIndex++] = num;
	m_nCurrentIndex %= m_nLength;
	return *this;
}

double TCircularBuffer::last()
{
	return data[(m_nCurrentIndex - 1 + m_nLength) % m_nLength];
}

int TCircularBuffer::length() const
{
	return this->m_nLength;
}

TCircularBuffer& TCircularBuffer::operator+=(double value)
{
	for (int i = 0; i < m_nLength; ++i)
		data[i] += value;
	return *this;
}

TCircularBuffer& TCircularBuffer::operator-=(double value)
{
	for (int i = 0; i < m_nLength; ++i)
		data[i] -= value;
	return *this;
}

TCircularBuffer& TCircularBuffer::operator*=(double value)
{
	for (int i = 0; i < m_nLength; ++i)
		data[i] *= value;
	return *this;
}

TCircularBuffer TCircularBuffer::operator*(double value)
{
	TCircularBuffer buffer(*this);
	for (int i = 0; i < m_nLength; ++i)
		buffer.data[i] *= value;
	return buffer;
}

TCircularBuffer& TCircularBuffer::operator/=(int value)
{
	for (int i = 0; i < m_nLength; ++i)
		data[i] /= value;
	return *this;
}

TCircularBuffer& TCircularBuffer::reverse()
{
	int idxA = m_nCurrentIndex,
		idxB = (m_nCurrentIndex - 1 + m_nLength) % m_nLength;
	while (idxA != idxB && idxA + 1 != idxB) {
		double t = data[idxA];
		data[idxA] = data[idxB];
		data[idxB] = t;
		idxA = idxA + 1 >= m_nLength ? 0 : idxA + 1;
		idxB = idxB - 1 < 0 ? m_nLength - 1 : idxB - 1;
		//std::cout << idxA << "\t" << idxB << std::endl;
	}
	return *this;
}

TCircularBuffer TCircularBuffer::toReverse()
{
	TCircularBuffer temp(m_nLength);
	temp.m_nCurrentIndex = 0;
	for (int i = 0; i < m_nLength; ++i)
		temp.data[i] = (*this)[m_nLength - 1 - i];
	return temp;
}

TCircularBuffer& TCircularBuffer::operator+=(const TCircularBuffer& buffer)
{
	if (buffer.m_nLength != m_nLength)
		throw("The lengthes of the buffers to be added must be the same.");
	int offset = (buffer.m_nCurrentIndex + m_nLength - m_nCurrentIndex) % m_nLength;
	double* dat = buffer.data;
	// although this implement is easy to comprehence, it's far beyond effective
	for (int i = 0; i < m_nLength - offset; ++i)
		data[i] += dat[i + offset];
	offset -= m_nLength;
	for (int i = -offset; i < m_nLength; ++i)
		data[i] += dat[i + offset];

	return *this;
}

TCircularBuffer TCircularBuffer::operator-()
{
	TCircularBuffer temp(*this);
	for (int i = 0; i < m_nLength; ++i)
		temp.data[i] = -data[i];
	return temp;
}

TCircularBuffer TCircularBuffer::slice(int start, int end) const
{
	start = start < 0 ? start + m_nLength : start;
	end = end < 0 ? end + m_nLength : end;
	int length = end - start + 1;
	length = length < 0 ? 0 : length;
	TCircularBuffer temp(length);
	for (int i = 0; i < length; ++i) {
		temp.data[i] = this->operator[](start + i);
	}
	return temp;
}

const double& TCircularBuffer::operator[](int index) const
{
	return data[(m_nCurrentIndex + index) % m_nLength];
}

TMatrix TCircularBuffer::toMatrix()
{
	TMatrix temp(m_nLength, 1);
	int index = 0;
	for (int i = m_nCurrentIndex; i < m_nLength; ++i)
		temp[index++][0] = data[i];
	for (int i = 0; i < m_nCurrentIndex; ++i)
		temp[index++][0] = data[i];
	return temp;
}

TMatrix TCircularBuffer::differ(bool bBackward)
{
	TMatrix temp(m_nLength - 1, 1);
	int sgn = bBackward ? 1 : -1;
	int index = 0;
	for (int i = m_nCurrentIndex; i < m_nLength - 1; ++i)
		temp[index++][0] = sgn * (data[i + 1] - data[i]);
	if (index < m_nLength - 1) temp[index++][0] = sgn * (data[0] - data[m_nLength - 1]);
	for (int i = 0; i < m_nCurrentIndex - 1; ++i)
		temp[index++][0] = sgn * (data[i + 1] - data[i]);
	return temp;
}

TCircularBuffer::const_iterator TCircularBuffer::begin()
{
	return const_iterator(*this);
}

TCircularBuffer::const_iterator TCircularBuffer::end()
{
	return const_iterator(*this, m_nLength);
}
