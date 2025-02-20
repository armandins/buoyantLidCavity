#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "src/fdOp.h"

// g++ -ggdb3 -O3 -std=c++23 -Wall -Wextra -pedantic -o sfvEnergyCavity45.out sfvEnergyCavity45.cc
// time ./sfvEnergyCavity45.out


void writeToVTK(const oneDVec& x, const oneDVec& y, 
                const twoDVec& u, 
                const twoDVec& v, 
                const twoDVec& Psi, 
                const twoDVec& Omega, 
                const twoDVec& Theta,
                const std::string& filename);

int main() {
    int   Nx              { 129           };              // Number of grid points in x-dir
    int   Ny              { 129           };              // Number of grid points in y-dir
    real  L               { 1.0           };              // Length of each wall
    real  UWall           { 1.0           };              // Upper-wall velocity
    real  dt              { 5E-05         };              // Timestep size
    real  maxe            { 1E-7          };              // Max-error
    real  Re              { 100           };              // Reynolds
    real  Pr              { 0.70          };              // Prandtl
    real  Gr              { 1E+3          };              // Grashof
    real  Pe              { Re*Pr         };              // Peclet
    real  Ra              { Gr / (Re * Re)};              // Rayleigh
    real  endTime         { 10.0          };             // End time for simulation
    real  time            { 0.0           };              // Current time

    twoDVec Omega  (Nx, oneDVec(Ny, 0.0)); // Vorticity
    twoDVec OmegaP (Nx, oneDVec(Ny, 0.0)); // Prev. Vorticity
    twoDVec Psi    (Nx, oneDVec(Ny, 0.0)); // Stream function
    twoDVec u      (Nx, oneDVec(Ny, 0.0)); // U velocity
    twoDVec v      (Nx, oneDVec(Ny, 0.0)); // V velocity
    twoDVec Theta  (Nx, oneDVec(Ny, 0.0)); // Theta 
    twoDVec ThetaP (Nx, oneDVec(Ny, 0.0)); // Prev. Theta

    // Grid setup
    real h = L / (Nx - 1);
    oneDVec x(Nx); 
    oneDVec y(Ny);
    for (int i = 0; i < Nx; ++i)    x[i] = i * h;
    for (int j = 0; j < Ny; ++j)    y[j] = j * h;

    // Time loop
    int iter = 0;
    while (time < endTime) {
        // Boundary Conditions for Psi - Omega
        for (int i = 0; i < Nx; ++i) {
            Omega[i][Ny - 1] = -2 * Psi[i][Ny - 2] / (h * h) - UWall * 2 / h;   // Top
            Omega[i][0] = -2 * Psi[i][1] / (h * h);                             // Bottom
        }
        for (int j = 0; j < Ny; ++j) {
            Omega[0][j] = -2 * Psi[1][j] / (h * h);                             // Left
            Omega[Nx - 1][j] = -2 * Psi[Nx - 2][j] / (h * h);                   // Right
        }

        // Boundary Conditions for Theta
        for (int j = 0; j < Ny; ++j) {
            Theta[0][j] = 1.0;                  // Left wall (Hot)
            Theta[Nx-1][j] = 0.0;               // Right wall (Cold)
        }
        for (int i = 0; i < Nx; ++i) {
            Theta[i][0] = Theta[i][1];          // Bottom wall (Zero gradient)
            Theta[i][Ny-1] = Theta[i][Ny-2];    // Top wall (Zero gradient)
        }

        // Solve Energy Equation (Explicit Euler)
        ThetaP = Theta;
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                real ddxTheta         = (ThetaP[i + 1][j] - ThetaP[i - 1][j]) / (2 * h);
                real ddyTheta         = (ThetaP[i][j + 1] - ThetaP[i][j - 1]) / (2 * h);
                real laplacianTheta   = (ThetaP[i + 1][j] + ThetaP[i - 1][j] + ThetaP[i][j + 1] + ThetaP[i][j - 1] - 4 * ThetaP[i][j]) / (h * h);
                real ddyPsi = (Psi[i][j + 1] - Psi[i][j - 1]) / (2 * h);
                real ddxPsi = (Psi[i + 1][j] - Psi[i - 1][j]) / (2 * h);

                Theta[i][j] = ThetaP[i][j] + (-ddyPsi * ddxTheta + ddxPsi * ddyTheta + (1.0 / Pe) * laplacianTheta) * dt;
            }
        }

        // Vorticity Transport Equation
        OmegaP = Omega;
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                real ddxTheta         = (ThetaP[i + 1][j] - ThetaP[i - 1][j]) / (2 * h);
                real ddyTheta         = (ThetaP[i][j + 1] - ThetaP[i][j - 1]) / (2 * h);
                real ddyPsi           = (Psi[i][j + 1] - Psi[i][j - 1]) / (2 * h);
                real ddxOmegaP        = (OmegaP[i + 1][j] - OmegaP[i - 1][j]) / (2 * h);
                real ddxPsi           = (Psi[i + 1][j] - Psi[i - 1][j]) / (2 * h);
                real ddyOmegaP        = (OmegaP[i][j + 1] - OmegaP[i][j - 1]) / (2 * h);
                real laplacianOmegaP  = (OmegaP[i + 1][j] + OmegaP[i - 1][j] + OmegaP[i][j + 1] + OmegaP[i][j - 1] - 4 * OmegaP[i][j]) / (h * h);
                real buoyancyTerm     = (sqrt(2))/2 * Ra * (ddxTheta - ddyTheta); 

                Omega[i][j] = OmegaP[i][j] - dt * (ddyPsi * ddxOmegaP - ddxPsi * ddyOmegaP) 
                            + dt * (1.0 / Re) * laplacianOmegaP 
                            + dt * buoyancyTerm;
            }
        }

        // Stream-function Laplacian
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                Psi[i][j] = (h * h * Omega[i][j] + Psi[i + 1][j] + Psi[i - 1][j] + Psi[i][j + 1] + Psi[i][j - 1]) / 4.0;
            }
        }

        // Create velocities from stream-function obtained above
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                u[i][j] = (Psi[i][j + 1] - Psi[i][j - 1]) / (2 * h);
                v[i][j] = (-Psi[i + 1][j] + Psi[i - 1][j]) / (2 * h);
            }
        }
        for (int i = 1; i < Nx - 1; ++i) {
            u[i][Ny - 1] = UWall; 
        }

        // Update time and iteration counter
        time += dt;
        iter++;

        // Convergence Check
        if (iter > 10) {
            real error = 0.0;
            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    error = std::max(error, std::abs(Omega[i][j] - OmegaP[i][j]));
                }
            }

            if (iter % 100 == 0)
                Info << "Iteration: " << iter << " | Time: " << time << " | Error: " << error << std::endl;
            if (error < maxe) {
                Info << "Converged after " << iter << " iterations." << std::endl;
                break;
            }
        }
    }

    // Final VTK output
    std::ostringstream vtkFilename;
    vtkFilename << "Gr" << std::fixed << std::setprecision(0) << Gr 
                << "Re" << std::fixed << std::setprecision(0) << Re 
                << "_final.vtk";
    writeToVTK(x, y, u, v, Psi, Omega, Theta, vtkFilename.str());

    return 0;
}

// ===============================================================================
// VTK Format Writer (*.vtk)                                                       
// ===============================================================================

void writeToVTK(const oneDVec& x, const oneDVec& y, 
                const twoDVec& u, 
                const twoDVec& v, 
                const twoDVec& Psi, 
                const twoDVec& Omega,
                const twoDVec& Theta,  
                const std::string&  filename) {
    std::ofstream vtkFile(filename);
    if (!vtkFile.is_open()) {
        cErr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    int Nx = x.size();
    int Ny = y.size();

    // VTK header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Lid-Driven Cavity Flow\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_GRID\n";
    vtkFile << "DIMENSIONS " << Nx << " " << Ny << " 1\n";
    vtkFile << "POINTS " << Nx * Ny << " double\n";

    // Rotation matrix components
    double cos45 = sqrt(2) / 2;
    double sin45 = sqrt(2) / 2;

    // Grid points with 45-degree rotation
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double x_rot = cos45 * x[i] - sin45 * y[j];
            double y_rot = sin45 * x[i] + cos45 * y[j];
            vtkFile << x_rot << " " << y_rot << " 0.0\n";
        }
    }

    // Write velocity data (vector field)
    vtkFile << "POINT_DATA " << Nx * Ny << "\n";
    vtkFile << "VECTORS velocity double\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            real u_rot = cos45 * u[i][j] - sin45 * v[i][j];
            real v_rot = sin45 * u[i][j] + cos45 * v[i][j];
            vtkFile << u_rot << " " << v_rot << " 0.0\n";
        }
    }

    // Write stream function (scalar field)
    vtkFile << "SCALARS stream_function double 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            vtkFile << Psi[i][j] << "\n";
        }
    }

    // Write u-velocity (scalar field)
    vtkFile << "SCALARS u_velocity double 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            vtkFile << u[i][j] << "\n";
        }
    }

    // Write v-velocity (scalar field)
    vtkFile << "SCALARS v_velocity double 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            vtkFile << v[i][j] << "\n";
        }
    }

    // Write vorticity (scalar field)
    vtkFile << "SCALARS vorticity double 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            vtkFile << Omega[i][j] << "\n";
        }
    }

    // Write theta (scalar field)
    vtkFile << "SCALARS theta double 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            vtkFile << Theta[i][j] << "\n";
        }
    }

    vtkFile.close();
    Info << "Results written to " << filename << std::endl;
}