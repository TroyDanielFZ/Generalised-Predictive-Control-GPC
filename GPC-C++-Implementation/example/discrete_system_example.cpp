#include <iostream>
#include <TSImulationSystem.h>
#include <TTrajectoryGenerator.h>
#include <vector>
#include <memory>

typedef std::vector<double> Vec;

void test_system(TSISO_Discrete& sys){
	{
		std::cout << "=== Step response ===" << std::endl;
		TTrajectory_Sin       traj (10.0, 0.1, 0.25);
		TTrajectory_Step control(1, 2);
		const double dt = 0.2;
		double t = 0.0;

		std::cout << "inputs: " << std::endl;
		while ((t += dt) < 10) {
			std::cout << control.getTrajectoryReference(t) << " ";
		}
		std::cout << std::endl << "outputs: " << std::endl;
		t = 0.0;
		while ((t += dt) < 10) {
			std::cout << sys.getOutput(control.getTrajectoryReference(t)) << " ";
		}
		std::cout << std::endl;
	}
	{
		std::cout << "=== Response of Rectangle Inputs ===" << std::endl;
		TTrajectory_Sin       traj (10.0, 0.1, 0.25);
		TTrajectory_Rectangle control(1, -1, 2);
		const double dt = 0.2;
		double t = 0.0;

		std::cout << "inputs: " << std::endl;
		while ((t += dt) < 10) {
			std::cout << control.getTrajectoryReference(t) << " ";
		}
		std::cout << std::endl << "outputs: " << std::endl;
		t = 0.0;
		while ((t += dt) < 10) {
			std::cout << sys.getOutput(control.getTrajectoryReference(t)) << " ";
		}
		std::cout << std::endl;
	}
}

int main() {
	// ----------------------------------------------------------------------------------
	// Test for TrajecotryGenerator
	// ----------------------------------------------------------------------------------

	// Contruct a system
	{
		std::cout << "================================================================================" << std::endl;
		std::cout << " test on a mininum-phase system" << std::endl;
		double num[]{ 1, 0.5};
		double den[]{ 1,  +1.5, 0.6 };
		Vec denominator(den, den + sizeof(den)/ sizeof(double));
		Vec numerator(num, num + sizeof(num)/ sizeof(double));

		TSISO_Discrete sys(denominator, numerator);
		test_system(sys);
		std::cout << "================================================================================" << std::endl;
	}
	{
		std::cout << "================================================================================" << std::endl;
		std::cout << " test on a non-mininum-phase system" << std::endl;
		double num[]{ 1, 0.5};
		double den[]{ 1,  -1.5, 0.6 };
		Vec denominator(den, den + sizeof(den)/ sizeof(double));
		Vec numerator(num, num + sizeof(num)/ sizeof(double));

		TSISO_Discrete sys(denominator, numerator);
		test_system(sys);
		std::cout << "================================================================================" << std::endl;
	}



	return 0;
}
