#pragma once

#include "VelocitySet.h"
#include "Subdomain.h"
#include "OutputFile.h"
#include "Boundary.h"
#include <vector>

class Simulation
{

public:
	int myrank_, num_procs_;
	SubDomain SubDomain_;
	Boundary Boundary_;
	double pressure_;
	double m_relaxation;
	std::vector<double> tau_w_x;
	std::vector<double> tau_w_y;
	std::vector<double> tau_w_z;
	std::vector<double> tau_w_magnitude;
	std::vector<double> InletVelocity_;
	Simulation(SubDomain& SubDomain, double tau, std::vector<double>& inVel, double pressure);
	~Simulation();
	void Collision();
	void Streaming();
	void ApplyInletEqulibriumBc();
	void ApplyInletZouHeBc();
	void ApplyOutletBc();
	void ApplyWallBc();
	double GetVelocityXGradient(std::shared_ptr<Node>& wall_node);
	double GetVelocityYGradient(std::shared_ptr<Node>& wall_node);
	double GetVelocityZGradient(std::shared_ptr<Node>& wall_node);
	void CalculateWallShearStress();
	void Communicate();
	void printData(const std::vector<double>& data, const std::string& label);
	void Iterate();
	void SetRelaxation(double tau);
	void WriteOutput(int iter);
};

