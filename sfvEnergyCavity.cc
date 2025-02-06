/*
=============================================================================
Arman Dindar Safa                                                      ++++++
2 / 6 / 2025                                                           ++++++
=============================================================================
This script computes Stream function, vorticity, velocity in x and y dir.
                     and Non-dimensionalized temperature for 
                     lid driven cavity problem using stream-function vor.
                     Continuity, momentum and energy equations are 
                     solved. Buoyant terms are taken into account.

Equations solved in latex format:
    I. Vorticity Transport Equation
    \frac{\partial \omega}{\partial t} + u \frac{\partial \omega}{\partial x} +
     v \frac{\partial \omega}{\partial y} = \frac{1}{Re} \nabla^2 \omega + 
     Ra \frac{\partial \Theta}{\partial x}

    II. Stream-function Laplacian

    III. Energy equation for Theta = T - Tc / Th - Tc 



Note: Navier-Stokes equations are non-dimensionalized, Which is why params
such as viscosity etc. are ommited from the equations.
Boundary conditions consist of: 
    I.      u/U = 1 Top wall
    II.     u/U = 0 Remaining walls
    III.    T - Tc / Th - Tc = Theta where Theta = 1 for left wall ( T = Th )
    IV.     Theta = 0 for right wall ( T = Tc )
    V.      Zero gradient for top and bottom walls.

In order to run the script, try: 

g++ -ggdb3 -O1 -std=c++23 -Wall -Wextra -pedantic -o sfvEnergyCavity.out sfvEnergyCavity.cc

Copyright: Use it as you wish, but in order to distribute the code or for 
modification purposes please keep my name in the script.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

void writeToVTK(const std::vector<double>& x, const std::vector<double>& y, 
                const std::vector<std::vector<double>>& u, 
                const std::vector<std::vector<double>>& v, 
                const std::vector<std::vector<double>>& Psi, 
                const std::vector<std::vector<double>>& Omega, 
                const std::vector<std::vector<double>>& Theta,
                const std::string& filename);

// Different test cases are for Gr = 1000, 10000 and 100000 to be computed

int main() {
    int     Nx              { 129   };              // Number of grid points in x-dir
    int     Ny              { 129   };              // Number of grid points in y-dir
    int     maxIt           { 500000};              // Maximum iter. allowed
    double  L               { 1.0   };              // Length of each wall
    double  UWall           { 1.0   };              // Upper-wall velocity
    double  dt              { 5e-05 };              // Timestep size
    double  maxe            { 1e-7  };              // Max-error
    double  Re              { 100   };              // Reynolds
    double  Pr              { 0.70  };              // Prandtl
    double  Gr              { 1000  };             // Grashof
    double  Pe              { Re*Pr };              // Peclet
    double  Ra              { Gr / (Re * Re)};      // Rayleigh


    std::vector<std::vector<double>> Omega  (Nx, std::vector<double>(Ny, 0.0)); // Vorticity
    std::vector<std::vector<double>> OmegaP (Nx, std::vector<double>(Ny, 0.0)); // Prev. Vorticity
    std::vector<std::vector<double>> Psi    (Nx, std::vector<double>(Ny, 0.0)); // Stream function
    std::vector<std::vector<double>> u      (Nx, std::vector<double>(Ny, 0.0)); // U velocity
    std::vector<std::vector<double>> v      (Nx, std::vector<double>(Ny, 0.0)); // V velocity
    std::vector<std::vector<double>> Theta  (Nx, std::vector<double>(Ny, 0.0)); // Theta 
    std::vector<std::vector<double>> ThetaP (Nx, std::vector<double>(Ny, 0.0)); // Prev. Theta


    //=============================================================================
    // Grid setup
    //-----------------------------------------------------------------------------
    double h = L / (Nx - 1);
    std::vector<double> x(Nx), y(Ny);
    for (int i = 0; i < Nx; ++i)    x[i] = i * h;
    for (int j = 0; j < Ny; ++j)    y[j] = j * h;


    // ======================================================
    // ======================================================
    // =====                SFV SOLUTION                =====
    // ======================================================
    // ======================================================
    
    for (int iter = 0; iter < maxIt; ++iter) {
        // ---------------------------------------------------------------------------
        // Boundary Conditions for Psi - Omega
        // ---------------------------------------------------------------------------
        for (int i = 0; i < Nx; ++i) {
            Omega[i][Ny - 1] = -2 * Psi[i][Ny - 2] / (h * h) - UWall * 2 / h; // Top
            Omega[i][0] = -2 * Psi[i][1] / (h * h); // Bottom
        }
        for (int j = 0; j < Ny; ++j) {
            Omega[0][j] = -2 * Psi[1][j] / (h * h); // Left
            Omega[Nx - 1][j] = -2 * Psi[Nx - 2][j] / (h * h); // Right
        }
        //-------------------------------------------------------------------------------
        // Boundary Conditions for Theta
        //-------------------------------------------------------------------------------
        for (int j = 0; j < Ny; ++j) {
            Theta[0][j] = 1.0;      // Left wall (Hot)
            Theta[Nx-1][j] = 0.0;   // Right wall (Cold)
        }
        for (int i = 0; i < Nx; ++i) {
            Theta[i][0] = Theta[i][1];          // Bottom wall (Zero gradient)
            Theta[i][Ny-1] = Theta[i][Ny-2];    // Top wall (Zero gradient)
        }

        //-------------------------------------------------------------------------------
        // Solve Energy Equation (Explicit Euler)
        //-------------------------------------------------------------------------------
        ThetaP = Theta;
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double ddxTheta         = (ThetaP[i + 1][j] - ThetaP[i - 1][j]) / (2 * h);
                double ddyTheta         = (ThetaP[i][j + 1] - ThetaP[i][j - 1]) / (2 * h);
                double laplacianTheta   = (ThetaP[i + 1][j] + ThetaP[i - 1][j] + ThetaP[i][j + 1] + ThetaP[i][j - 1] - 4 * ThetaP[i][j]) / (h * h);
                double ddyPsi = (Psi[i][j + 1] - Psi[i][j - 1]) / (2 * h);
                double ddxPsi = (Psi[i + 1][j] - Psi[i - 1][j]) / (2 * h);
                //Theta[i][j]             = ThetaP[i][j] - dt * (u[i][j] * ddxTheta + v[i][j] * ddyTheta) 
                //                        + dt * (1.0 / Pe) * laplacianTheta;
                Theta[i][j] = ThetaP[i][j] + (-ddyPsi * ddxTheta + ddxPsi * ddyTheta + (1.0 / Pe) * laplacianTheta) * dt;
            }
        }


        //----------------------------------------------------------------------------
        // Vorticity Transport Equation
        //----------------------------------------------------------------------------
        OmegaP = Omega;
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double ddxTheta = (ThetaP[i + 1][j] - ThetaP[i - 1][j]) / (2 * h);
                double ddyPsi = (Psi[i][j + 1] - Psi[i][j - 1]) / (2 * h);
                double ddxOmegaP = (OmegaP[i + 1][j] - OmegaP[i - 1][j]) / (2 * h);
                double ddxPsi = (Psi[i + 1][j] - Psi[i - 1][j]) / (2 * h);
                double ddyOmegaP = (OmegaP[i][j + 1] - OmegaP[i][j - 1]) / (2 * h);
                double laplacianOmegaP = (OmegaP[i + 1][j] + OmegaP[i - 1][j] + OmegaP[i][j + 1] + OmegaP[i][j - 1] - 4 * OmegaP[i][j]) / (h * h);
                double buoyancyTerm = Ra * ddxTheta; 

                Omega[i][j] = OmegaP[i][j] - dt * (ddyPsi * ddxOmegaP - ddxPsi * ddyOmegaP) 
                            + dt * (1.0 / Re) * laplacianOmegaP 
                            + dt * buoyancyTerm;
            }
        }

        //----------------------------------------------------------------------------
        // Stream-function Laplacian
        //----------------------------------------------------------------------------
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                Psi[i][j] = (h * h * Omega[i][j] + Psi[i + 1][j] + Psi[i - 1][j] + Psi[i][j + 1] + Psi[i][j - 1]) / 4.0;
            }
        }

        //----------------------------------------------------------------------------
        // Convergence Check
        //----------------------------------------------------------------------------
        if (iter > 10) {
            double error = 0.0;
            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    error = std::max(error, std::abs(Omega[i][j] - OmegaP[i][j]));
                }
            }

                // Print error every 100 iterations
            if (iter % 100 == 0)
                std::cout << "Iteration: " << iter << " | Error: " << error << std::endl;
            if (error < maxe) {
                std::cout << "Converged after " << iter << " iterations." << std::endl;
                break;
            }
        }
    }
    //----------------------------------------------------------------------------
    // Create velocities from stream-function obtained above
    //----------------------------------------------------------------------------
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            u[i][j] = (Psi[i][j + 1] - Psi[i][j - 1]) / (2 * h);
            v[i][j] = (-Psi[i + 1][j] + Psi[i - 1][j]) / (2 * h);
        }
    }
    for (int i = 1; i < Nx - 1; ++i) {
        u[i][Ny - 1] = UWall; 
    }

    //----------------------------------------------------------------------------
    // Call VTK Writer
    //----------------------------------------------------------------------------
    writeToVTK(x, y, u, v, Psi, Omega, Theta, "Gr1000cavityResults.vtk");
    return 0;
}

// ===============================================================================
// VTK Format Writer (*.vtk)                                                       
// ===============================================================================

void writeToVTK(const std::vector<double>& x, const std::vector<double>& y, 
                const std::vector<std::vector<double>>& u, 
                const std::vector<std::vector<double>>& v, 
                const std::vector<std::vector<double>>& Psi, 
                const std::vector<std::vector<double>>& Omega,
                const std::vector<std::vector<double>>& Theta,  
                const std::string& filename) {
    std::ofstream vtkFile(filename);
    if (!vtkFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
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

    // Grid points
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            vtkFile << x[i] << " " << y[j] << " 0.0\n";
        }
    }

    // Write velocity data (vector field)
    vtkFile << "POINT_DATA " << Nx * Ny << "\n";
    vtkFile << "VECTORS velocity double\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            vtkFile << u[i][j] << " " << v[i][j] << " 0.0\n";
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
    std::cout << "Results written to " << filename << std::endl;
}

// Might add a TECPLOT writer aswell as a python post-processor aswell later on.