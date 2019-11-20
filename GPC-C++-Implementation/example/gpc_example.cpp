#include <iostream>
#include <iomanip>
#include "GPC.h"
#include "TSImulationSystem.h"
#include "TTrajectoryGenerator.h"
#include <vector>
#include <memory>
#include <fstream>

// Pre-define and constants
typedef std::vector<double> Vec;
using std::setw;
const double dt = 0.2;

void gpc_control(){

	// ----------------------------------------------------------------------------------
	// Test for TrajecotryGenerator
	// ----------------------------------------------------------------------------------
	std::freopen("log.txt", "w", stdout);
	TTrajectory_Sin traj (1.0, 0.2, pi /4 );
	TTrajectory_Rectangle control(1, -1, 2);

	// Construct a system
	double num[]{ 1, 0.5};
	double den[]{ 1,  +1.5, 0.6 };
	Vec denominator(den, den + sizeof(den)/ sizeof(double));
	Vec numerator(num, num + sizeof(num)/ sizeof(double));
	TSISO_Discrete sys(denominator, numerator);

	// All the parameters you can tune here
	int N = 5, Nu = 3, na = 3, nb = 2;
	double alpha = 0.3, lambda = 0.4, rho = 0.99;
	double beta = 1.5;
	// construct a gpc controller
	std::unique_ptr<GPC> gpc = std::unique_ptr<GPC>(new BetaGPC(na, nb, N, Nu, rho, alpha, lambda, beta));

	double t = 0.0;
	double y = 0;

	// try to load the parameters, but no check whether the format is correct of not
	gpc->loadParameters("param.txt"); 
	while (t < 800) {
		double yr = traj.getTrajectoryReference(t);
		double u = gpc->getControl(yr, y);
		y = sys.getOutput(u);
		std::cout << setw(8)<< t << ", " << setw(8)<< yr << ", " << setw(8)<< u << ", "<< setw(8)<< y << ";" << std::endl;
		t += dt;
	}
	gpc->saveParameters("param.txt"); // and save the parameters for further analysis/usage
}

int main() {
	gpc_control();
	return 0;
}
