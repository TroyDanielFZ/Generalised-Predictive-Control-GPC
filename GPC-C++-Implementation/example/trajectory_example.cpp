#include <iostream>
#include <TTrajectoryGenerator.h>
#include <memory>


int main(){
	// ----------------------------------------------------------------------------------
	// Test for TrajecotryGenerator
	// ----------------------------------------------------------------------------------
	{
		TTrajectoryGenerator traj = TTrajectory_Sin(10.0, 0.1, 0.25);
		const double dt = 0.2;
		for (int i = 0; i < 100; ++i) {
			std::cout << traj.getTrajectoryReference(dt * i) << " ";
		}
		std::cout << std::endl 
			<< "------------------------------------------------------------" << std::endl;
	}
	{
		const TTrajectoryGenerator& traj = TTrajectory_Sin(10.0, 0.1, 0.25);
		const double dt = 0.2;
		for (int i = 0; i < 100; ++i) {
			std::cout << traj.getTrajectoryReference(dt * i) << " ";
		}
		std::cout << std::endl 
			<< "------------------------------------------------------------" << std::endl;
	}
	{ // By using pointer, it works properly
		// And, by using pointer, the base class can be an virutal 
		// interface, e.g. pure base class
		TTrajectoryGenerator * traj =  new TTrajectory_Sin(10.0, 0.1, 0.25);
		const double dt = 0.2;
		for (int i = 0; i < 100; ++i) {
			std::cout << traj->getTrajectoryReference(dt * i) << " ";
		}
		delete traj;
		std::cout << std::endl 
			<< "------------------------------------------------------------" << std::endl;
	}
	{
		std::unique_ptr<TTrajectoryGenerator> traj(new TTrajectory_Sin(10.0, 0.1, 0.25));
		const double dt = 0.2;
		for (int i = 0; i < 100; ++i) {
			std::cout << traj->getTrajectoryReference(dt * i) << " ";
		}
		std::cout << std::endl 
			<< "------------------------------------------------------------" << std::endl;
	}
}
