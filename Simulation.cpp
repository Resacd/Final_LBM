#include "Simulation.h"

Simulation::Simulation(SubDomain& SubDomain,double tau,std::vector<double> &inVel,double pressure)
	: SubDomain_(SubDomain)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs_);
	SetRelaxation(tau);
    InletVelocity_ = inVel;
    pressure_ = pressure;
    //std::cout << "Processor " << myrank_ << ": Nx start = " << SubDomain_.start_x_ << ", Nx end = " << SubDomain_.end_x_ << std::endl;
    //std::cout << "Processor " << myrank_ << ": Ny start = " << SubDomain_.start_y_ << ", Ny end = " << SubDomain_.end_y_ << std::endl;
    //std::cout << "Processor " << myrank_ << ": Nz start = " << SubDomain_.start_z_ << ", Nz end = " << SubDomain_.end_z_ << std::endl;
    Iterate();
}

Simulation::~Simulation()
{
}


void Simulation::Collision()
{
    for (int k = 1; k < SubDomain_.my_Nz_-1; k++) {
        for (int j = 1; j < SubDomain_.my_Ny_-1; j++) {
            for (int i = 1; i < SubDomain_.my_Nx_-1; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {

                    if (SubDomain_.lattice_[i][j][k]!=nullptr) {

                        SubDomain_.lattice_[i][j][k]->m_newDistributions[dir] =
                            SubDomain_.lattice_[i][j][k]->m_distributions[dir] -
                            (1.0 / m_relaxation) * (SubDomain_.lattice_[i][j][k]->m_distributions[dir]-
                                SubDomain_.lattice_[i][j][k]->Equilibrium(SubDomain_.GetVelSet(), dir));

                    }
                }
            }
        }
    }
}


void Simulation::Streaming()
{
    for (int k = 0; k < SubDomain_.my_Nz_; k++) {
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {

                    if (SubDomain_.lattice_[i][j][k]!=nullptr && SubDomain_.lattice_[i][j][k]->m_neighbours[dir] != nullptr) {
                        SubDomain_.lattice_[i][j][k]->Stream(dir);
                    }
                }
            }
        }
    }
}


void Simulation::printData(const std::vector<double>& data, const std::string& label) {
    std::cout << "Processor rank: " << myrank_ << std::endl;
    std::cout << label << " data:"<<std::endl;
    for (const auto& value : data) {
        std::cout << value<<std::endl;
    }
    std::cout << std::endl;
}

void Simulation::Communicate() {

    int tag_xp = 0;
    int tag_xm = 1;
    int tag_yp = 2;
    int tag_ym = 3;
    int tag_zp = 4;
    int tag_zm = 5;

    MPI_Datatype Distribution;
    MPI_Type_contiguous(_nLatNodes, MPI_DOUBLE, &Distribution);
    MPI_Type_commit(&Distribution);

    MPI_Status status;

    if (SubDomain_.my_right_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    if (SubDomain_.lattice_[SubDomain_.my_Nx_ - 2][j][k] == nullptr) {
                        send_data.push_back(0.0);
                    }
                    else {

                        send_data.push_back(SubDomain_.lattice_[SubDomain_.my_Nx_ - 2][j][k]->m_distributions[dir]);
                    }
                    
                }
            }
        }
        //std::cout << std::endl;
       // printf("Process %d: Sending data to the right neighbor\n", myrank_);


        std::vector<double> recv_data(send_data.size());

        //printData(send_data, "Send (right)");

        MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_right_, tag_xp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_right_, tag_xm,
            MPI_COMM_WORLD, &status);
        //std::cout << std::endl;
       // printf("Process %d: Received data from the right neighbor\n", myrank_);

        //printData(recv_data, "Receive (right)");

        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {

                    SubDomain_.lattice_[SubDomain_.my_Nx_ - 1][j][k]->m_distributions[dir] = recv_data[index++];
                }
            }
        }
    }


    if (SubDomain_.my_left_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {

                    if (SubDomain_.lattice_[1][j][k] == nullptr) {
                        send_data.push_back(0.0);
                    }
                    else {

                        send_data.push_back(SubDomain_.lattice_[1][j][k]->m_distributions[dir]);

                    }
                }
            }
        }
        //std::cout << std::endl;
       // printf("Process %d: Sending data to the left neighbor\n", myrank_);
        std::vector<double> recv_data(send_data.size());
        //printData(send_data, "Send (left)");
        MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_left_, tag_xm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_left_, tag_xp,
            MPI_COMM_WORLD, &status);
        //printData(recv_data, "Receive (left)");
        //std::cout << std::endl;
       // printf("Process %d: Received data from the left neighbor\n", myrank_);
        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.lattice_[0][j][k]->m_distributions[dir] = recv_data[index++];
                }
            }
        }
    }

    if (SubDomain_.my_front_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {

                    if (SubDomain_.lattice_[i][1][k] == nullptr) {

                        send_data.push_back(0.0);
                    }
                    else {

                        send_data.push_back(SubDomain_.lattice_[i][1][k]->m_distributions[dir]);
                    }
                }
            }
        }
       // printf("Process %d: Sending data to the front neighbor\n", myrank_);
        std::vector<double> recv_data(send_data.size(), 0.0);
        //printData(send_data, "Send (front)");
        MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_front_, tag_ym,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_front_, tag_yp,
            MPI_COMM_WORLD, &status);
       // printf("Process %d: Received data from the front neighbor\n", myrank_);
        //printData(recv_data, "Receive (front)");
        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.lattice_[i][0][k]->m_distributions[dir] = recv_data[index++];
                }
            }
        }
    }

    if (SubDomain_.my_back_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {

                    if (SubDomain_.lattice_[i][SubDomain_.my_Ny_ - 2][k] == nullptr) {
                        send_data.push_back(0.0);
                    }
                    else {
                        send_data.push_back(SubDomain_.lattice_[i][SubDomain_.my_Ny_ - 2][k]->m_distributions[dir]);
                    }
                }
            }
        }
       // printf("Process %d: Sending data to the back neighbor\n", myrank_);
        std::vector<double> recv_data(send_data.size(), 0.0);
        //printData(send_data, "Send (back)");
        MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_back_, tag_yp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_back_, tag_ym,
            MPI_COMM_WORLD, &status);

       // printf("Process %d: Received data from the back neighbor\n", myrank_);
        //printData(recv_data, "Receive (back)");
        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.lattice_[i][SubDomain_.my_Ny_ - 1][k]->m_distributions[dir] = recv_data[index++];
                }
            }
        }
    }


    if (SubDomain_.my_up_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {

                    if (SubDomain_.lattice_[i][j][SubDomain_.my_Nz_ - 2] == nullptr) {

                        send_data.push_back(0.0);
                    }
                    else {

                        send_data.push_back(SubDomain_.lattice_[i][j][SubDomain_.my_Nz_ - 2]->m_distributions[dir]);
                    }
                }
            }
        }
       // printf("Process %d: Sending data to the up neighbor\n", myrank_);
        std::vector<double> recv_data(send_data.size(), 0.0);
        //printData(send_data, "Send (up)");
        MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_up_, tag_zp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_up_, tag_zm,
            MPI_COMM_WORLD, &status);
       // printf("Process %d: Received data from the up neighbor\n", myrank_);
        //printData(recv_data, "Receive (up)");
        int index = 0;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.lattice_[i][j][SubDomain_.my_Nz_ - 1]->m_distributions[dir] = recv_data[index++];
                }
            }
        }
    }


    if (SubDomain_.my_down_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {

                    if (SubDomain_.lattice_[i][j][1] == nullptr) {

                        send_data.push_back(0.0);
                    }

                    else {

                        send_data.push_back(SubDomain_.lattice_[i][j][1]->m_distributions[dir]);
                    }
                }
            }
        }
       // printf("Process %d: Sending data to the down neighbor\n", myrank_);
        std::vector<double> recv_data(send_data.size(), 0.0);
        //printData(send_data, "Send (down)");
        MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_down_, tag_zm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_down_, tag_zp,
            MPI_COMM_WORLD, &status);
       // printf("Process %d: Received data from the down neighbor\n", myrank_);
        //printData(recv_data, "Receive (down)");
        int index = 0;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.lattice_[i][j][0]->m_distributions[dir] = recv_data[index++];
                }
            }
        }
    }

}


void Simulation::Iterate()
{
    for (int iter = 0; iter < 50; iter++) {
        Collision();
        Streaming();
        ApplyInletEqulibriumBc();
        ApplyOutletBc();
        ApplyWallBc();
        Communicate();
        WriteOutput(iter);
        if (myrank_ == 0) {
            std::cout << std::endl;
            std::cout << "Iteration " << " " << iter << " " << "has been performed successfully " << std::endl;
        }
    }
}


void Simulation::ApplyInletEqulibriumBc() {

    for (int lat = 0; lat < SubDomain_.InletBoundaryNode_.size(); lat++) {
        for (int dir = 0; dir < _nLatNodes; dir++) {
            Boundary_.Apply_Equilibrium_Bc(SubDomain_.velSet_, *(SubDomain_.InletBoundaryNode_[lat])
                , InletVelocity_, dir);
        }
    }
}

void Simulation::ApplyInletZouHeBc() {

    for (int lat = 0; lat < SubDomain_.InletBoundaryNode_.size(); lat++) {

        for (int dir = 0; dir < _nLatNodes; dir++) {
            Boundary_.Apply_ZouHe_Bc(SubDomain_.velSet_, *(SubDomain_.InletBoundaryNode_[lat])
                , InletVelocity_, dir);
        }
    }
}

void Simulation::ApplyOutletBc() {

    for (int lat = 0; lat < SubDomain_.OutletBoundaryNode_.size(); lat++) {
        for (int dir = 1; dir < _nLatNodes; dir++)
        {
            Boundary_.ApplyConstPressure(SubDomain_.velSet_, *(SubDomain_.OutletBoundaryNode_[lat]), dir);
        }
    }
}

void Simulation::ApplyWallBc() {

    for (int lat = 0; lat < SubDomain_.WallBoundaryNode_.size(); lat++) {
        for (int dir = 1; dir < _nLatNodes; dir++) {

            int opp_dir = SubDomain_.GetVelSet()->OppositeDirection(dir);
            if (SubDomain_.WallBoundaryNode_[lat]->m_neighbours[dir] == nullptr) {

                Boundary_.Apply_BounceBack(*(SubDomain_.WallBoundaryNode_[lat]), dir, opp_dir);

            }
        }
    }

}

double Simulation::GetVelocityXGradient(std::shared_ptr<Node>& wall_node) {
    double gradient = 0.0;

    for (int dir = 0; dir < _nLatNodes; dir++) {
        if (wall_node->m_neighbours[dir] != nullptr) {
            double delta_x = SubDomain_.Geometry_.spacing;
            double delta_vx = wall_node->m_neighbours[dir]->Velocity(SubDomain_.velSet_)[0]
                - wall_node->Velocity(SubDomain_.velSet_)[0];
            gradient += delta_vx / delta_x;
        }
    }
    return gradient;
}

double Simulation::GetVelocityYGradient(std::shared_ptr<Node>& wall_node) {
    double gradient = 0.0;

    for (int dir = 0; dir < _nLatNodes; dir++) {
        if (wall_node->m_neighbours[dir] != nullptr) {
            double delta_y = SubDomain_.Geometry_.spacing;
            double delta_vy = wall_node->m_neighbours[dir]->Velocity(SubDomain_.velSet_)[1]
                - wall_node->Velocity(SubDomain_.velSet_)[1];
            gradient += delta_vy / delta_y;
        }
    }
    return gradient;
}


double Simulation::GetVelocityZGradient(std::shared_ptr<Node>& wall_node) {
    double gradient = 0.0;

    for (int dir = 0; dir < _nLatNodes; dir++) {
        if (wall_node->m_neighbours[dir] != nullptr) {
            double delta_z = SubDomain_.Geometry_.spacing;
            double delta_vz = wall_node->m_neighbours[dir]->Velocity(SubDomain_.velSet_)[2]
                - wall_node->Velocity(SubDomain_.velSet_)[2];
            gradient += delta_vz / delta_z;
        }
    }
    return gradient;
}


void Simulation::CalculateWallShearStress() {
    double mu = 1e-6;

    tau_w_x.resize(SubDomain_.WallBoundaryNode_.size());
    tau_w_y.resize(SubDomain_.WallBoundaryNode_.size());
    tau_w_z.resize(SubDomain_.WallBoundaryNode_.size());
    tau_w_magnitude.resize(SubDomain_.WallBoundaryNode_.size());

    for (int lat = 0; lat < SubDomain_.WallBoundaryNode_.size(); lat++) {
        double dudx = GetVelocityXGradient(SubDomain_.WallBoundaryNode_[lat]);
        double dudy = GetVelocityYGradient(SubDomain_.WallBoundaryNode_[lat]);
        double dudz = GetVelocityZGradient(SubDomain_.WallBoundaryNode_[lat]);

        tau_w_x[lat] = mu * dudx;
        tau_w_y[lat] = mu * dudy;
        tau_w_z[lat] = mu * dudz;

        tau_w_magnitude[lat] = sqrt(tau_w_x[lat] * tau_w_x[lat]
            + tau_w_y[lat] * tau_w_y[lat] + tau_w_z[lat] * tau_w_z[lat]);
    }
}



void Simulation::SetRelaxation(double tau)
{
	m_relaxation = tau;
}

void Simulation::WriteOutput(int iter)
{
	OutputFile file(iter);
	file.WriteCSV(SubDomain_.lattice_, SubDomain_.velSet_,
        SubDomain_.start_x_,SubDomain_.start_y_,SubDomain_.start_z_,
        SubDomain_.end_x_,SubDomain_.end_y_,SubDomain_.end_z_);

}