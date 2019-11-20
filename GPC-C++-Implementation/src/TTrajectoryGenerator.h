#pragma once
class TTrajectoryGenerator {
private:
public:
	TTrajectoryGenerator();
	~TTrajectoryGenerator();
	virtual double getTrajectoryReference(double t) const;
};

class TTrajectory_Sin : public TTrajectoryGenerator {
	double _frequency;
	double _phase_lag;
	double _amplitude;
public:
	TTrajectory_Sin(double amplitdue, double frequency, double phase_lag);
	double getTrajectoryReference(double t) const override;
};

class TTrajectory_Step : public TTrajectoryGenerator {
	double _amplitude;
	double _step_time;
public:
	TTrajectory_Step(double amplitdue, double step_time);
	double getTrajectoryReference(double t) const override;
};
class TTrajectory_Rectangle : public TTrajectoryGenerator {
	double _amplitude_pos;
	double _amplitude_neg;
	double _half_period;
public:
	TTrajectory_Rectangle(double amplitdue_pos, double amplitdue_neg, double period);
	double getTrajectoryReference(double t) const override;
};