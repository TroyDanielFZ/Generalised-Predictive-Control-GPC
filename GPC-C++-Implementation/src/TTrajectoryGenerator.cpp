#include "TTrajectoryGenerator.h"
#include <cmath>
#include "TMath.h"

TTrajectoryGenerator::TTrajectoryGenerator() { }

TTrajectoryGenerator::~TTrajectoryGenerator() { }

double TTrajectoryGenerator::getTrajectoryReference(double t)  const{
	return 0.0;
}

TTrajectory_Sin::TTrajectory_Sin(double amplitdue, double frequency, double phase_lag)
	: _amplitude(amplitdue)
	, _frequency(frequency)
	, _phase_lag(phase_lag) { }

double TTrajectory_Sin::getTrajectoryReference(double t)  const{
	//return _amplitude* sin(2*pi * _frequency * t + _phase_lag);
	return _amplitude* sin( _frequency * t + _phase_lag);
}

TTrajectory_Step::TTrajectory_Step(double amplitdue, double step_time)
	:_amplitude(amplitdue)
	, _step_time(step_time){ }

double TTrajectory_Step::getTrajectoryReference(double t) const {
	return (t>= _step_time) ? _amplitude: 0.0;
}

TTrajectory_Rectangle::TTrajectory_Rectangle(double amplitdue_pos, double amplitdue_neg, double period)
	: _amplitude_pos(amplitdue_pos)
	, _amplitude_neg(amplitdue_neg)
	, _half_period(period/ 2) { }

double TTrajectory_Rectangle::getTrajectoryReference(double t)  const{
	return int((t - _half_period * .0000001)/_half_period) % 2 ? _amplitude_pos : _amplitude_neg;
}
