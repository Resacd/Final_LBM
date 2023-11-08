#pragma once
#include <vector>
#include <iostream>
#include "Struct.h"


int gridx = 25;
double U = 0.1;
double Re = 100;
double niu = U * gridx / Re;
double tau = 3.0 * niu + 0.5;
double spacing = 0.1;
int mdepth = 1;
std::vector<double> inVel{ 1.0,0.0,0.0 };
std::vector<double> InletCtr{ 36.441,-208.497,167.44 };
double Inletradius = 1.126;
std::vector<double> OutletCtr{ 14.2978,-211.996,-172.103 };
double Outletradius = 1.5235;
double pressure = 1.0;