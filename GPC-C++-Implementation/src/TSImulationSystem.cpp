#include "TSImulationSystem.h"

TSISO::TSISO() { }

TSISO::~TSISO() { }


TSISO_Discrete::TSISO_Discrete(std::vector<double> denominator, std::vector<double> numerator)
	:TSISO() 
	, _denominator(denominator.size() - 1)
	, _numerator(numerator)
	, _y_history(_denominator.size())
	, _u_history(_numerator.size()){
	// Confirm the ratio
	//auto it = denominator.rbegin();
	auto it = denominator.begin();
	double ratio = *it;
	while (abs(ratio) < 1e-6) (ratio = *(++it));
	_ASSERT(ratio >= 1e-6);
	it =  ++ denominator.begin();
	auto it_to = _denominator.begin();
	while (it < denominator.end()) { *it_to ++ = *it ++ / ratio; }
	it = _numerator.begin();
	while (it < _numerator.end()) { *it ++ /= ratio; }
}

TSISO_Discrete::~TSISO_Discrete() { }

double TSISO_Discrete::getOutput(double input) {
	_u_history.push(input);
	double y = 0.0;

	auto len= _y_history.length();
	for (int i = 0; i < len; ++i)
		y -= _denominator[len - 1UL - i] * _y_history[i];
	len = _u_history.length() - 1;
	for (auto i = 0UL; i <= len; ++i)
		y += _numerator[len - i] * _u_history[i];
	_y_history.push(y);
	return y;
}


std::ostream& operator<<(std::ostream& os, const TSISO_Discrete& sys) {
	std::cout << std::endl << "This system is :" << std::endl;
	int len = sys._denominator.size();
	std::cout << sys._numerator[0] << "z(" << ((int) sys._numerator.size() -len - 1) << ")";
	for (int i = sys._numerator.size() - 1; i > 0; --i) {
		std::cout << " + "<< sys._numerator[i] << "z(" << (i - 1- len) << ") ";
	}
	std::cout << std::endl  << "1";
	for (int i = 0; i < sys._denominator.size(); ++i) {
		std::cout << " + "<< sys._denominator[i] << "z(-" << i + 1  << ") ";
	}
	return os;
}
