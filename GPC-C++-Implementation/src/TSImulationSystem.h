#pragma once
#include <iostream>
#include <vector>
#include "TMath.h"
class TSISO {
public:
	TSISO();
	~TSISO();
	virtual double getOutput(double input) = 0;
};
// ----------------------------------------------------------------------------------
// The system is formulated as:
//      B(z)      b0 + b0z^-1 + ...
// y = ----- u = ------------------------ u, e.g. A(z)y = B(z) u
//      A(Z)      a0 + a0z^-1 + ...
// A(z) is in form [a0, a1, a2, ...]
// B(z) is in form [b0, b1, b2, ...]
// ----------------------------------------------------------------------------------
class TSISO_Discrete : public TSISO {
private:
	std::vector<double> _denominator;
	std::vector<double> _numerator;
	TCircularBuffer _y_history;
	TCircularBuffer _u_history;
	
public:
	TSISO_Discrete(std::vector<double> denominator, std::vector<double> numerator);
	~TSISO_Discrete();
	double getOutput(double input);
	friend std::ostream& operator<< (std::ostream& os, const TSISO_Discrete&);
};

