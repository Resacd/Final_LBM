#include <iostream>
#include "LatticeModel.h"
#include "Simulation.h"
#include "OutputFile.h"
#include <mpi.h>
#include "Subdomain.h"
#include "Subdomain.h"
#include "Input.h"
#include "ComplexGeometry.h"
#include "disclaimer.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <iomanip>
#include <omp.h>
#include "Utility.h"
#include "ImageProcessor.h"
#include "ConfigReader.h"

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    int myrank, num_procs;
    Utility Utility;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    if (myrank == 0) {
        Print_disclaimer();
        //ConfigReader config("Input.txt");
        //std::map<std::string, double> configData;
        //std::string imageFile = "inlet.png";
        //ImageProcessor ImageProcessor(imageFile);
        ////ImageProcessor.displayImage();

        //if (config.readConfigData(configData)) {
        //    double minVFR = configData["minVFR"];
        //    double maxVFR = configData["maxVFR"];
        //    double minTime = configData["minTime"];
        //    double maxTime = configData["maxTime"];

        //    std::vector<double> vfrData;
        //    std::vector<double> timeData;

        //    ImageProcessor.extractDataFromImage(imageFile, minVFR, maxVFR, minTime, maxTime, vfrData, timeData);
        //    std::string velocityEquation = ImageProcessor.fitVelocityEquation(vfrData, timeData);
        //    std::cout << "Obtained Velocity Equation: " << velocityEquation << std::endl;
        //}



        omp_set_num_threads(omp_get_max_threads());
        int numThreads = omp_get_max_threads();
        omp_set_num_threads(numThreads);
        bool stopChrono = false;
        std::thread chronoThread([&]() {
            while (!stopChrono) {
                auto now = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = now - start;

                std::cout << "\rElapsed time: " << std::fixed << std::setprecision(2)
                    << std::setw(8) << elapsed.count() << " seconds" << std::flush;

                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
            });
        chronoThread.detach();
    }

    std::shared_ptr<LatticeModel> velSet = std::make_shared<LatticeModel>();

    Geometry Geometry("CylinderASCIISTL.stl", InletCtr, Inletradius, OutletCtr, Outletradius, spacing, mdepth);
    SubDomain SubDomain(velSet, Geometry);
    Simulation(SubDomain, tau, inVel, pressure);

    MPI_Finalize();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    if (myrank == 0) {
        std::cout << "\rElapsed time: " << std::fixed << std::setprecision(2)
            << std::setw(8) << elapsed.count() << " seconds." << std::endl;
    }

    return 0;
}
