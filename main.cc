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
     Ri \frac{\partial \Theta}{\partial x}

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

or simply type:
chmod +x run.sh
./run.sh 1
Where 1 is the optimization level. Choose a number between 1, 2 and 3. 

Copyright: Use it as you wish, but in order to distribute the code or for 
modification purposes please keep my name in the script.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "src/fdOp.h"
#include "src/vtkWriter.h"
#include "src/constParams.h"

// Different test cases are for Gr = 1000, 10000 and 100000 to be computed

int main() {
    
    twoDVec Omega  (Nx, oneDVec(Ny, 0.0)); // Vorticity
    twoDVec OmegaP (Nx, oneDVec(Ny, 0.0)); // Prev. Vorticity
    twoDVec Psi    (Nx, oneDVec(Ny, 0.0)); // Stream function
    twoDVec u      (Nx, oneDVec(Ny, 0.0)); // U velocity
    twoDVec v      (Nx, oneDVec(Ny, 0.0)); // V velocity
    twoDVec Theta  (Nx, oneDVec(Ny, 0.0)); // Theta 
    twoDVec ThetaP (Nx, oneDVec(Ny, 0.0)); // Prev. Theta


    //=============================================================================
    // Grid setup
    //-----------------------------------------------------------------------------
    real h = L / (Nx - 1);
    oneDVec x(Nx), y(Ny);
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
            Omega[i][Ny - 1] = -2 * Psi[i][Ny - 2] / (h * h) - UWall * 2 / h;   // Top
            Omega[i][0] = -2 * Psi[i][1] / (h * h);                             // Bottom
        }
        for (int j = 0; j < Ny; ++j) {
            Omega[0][j] = -2 * Psi[1][j] / (h * h);                             // Left
            Omega[Nx - 1][j] = -2 * Psi[Nx - 2][j] / (h * h);                   // Right
        }
        //-------------------------------------------------------------------------------
        // Boundary Conditions for Theta
        //-------------------------------------------------------------------------------
        for (int j = 0; j < Ny; ++j) {
            Theta[0][j] = 1.0;                  // Left wall (Hot)
            Theta[Nx-1][j] = 0.0;               // Right wall (Cold)
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
                real ddxTheta         = (ThetaP[i + 1][j] - ThetaP[i - 1][j]) / (2 * h);
                real ddyTheta         = (ThetaP[i][j + 1] - ThetaP[i][j - 1]) / (2 * h);
                real laplacianTheta   = (ThetaP[i + 1][j] + ThetaP[i - 1][j] + ThetaP[i][j + 1] + ThetaP[i][j - 1] - 4 * ThetaP[i][j]) / (h * h);
                real ddyPsi = (Psi[i][j + 1] - Psi[i][j - 1]) / (2 * h);
                real ddxPsi = (Psi[i + 1][j] - Psi[i - 1][j]) / (2 * h);

                Theta[i][j] = ThetaP[i][j] + (-ddyPsi * ddxTheta + ddxPsi * ddyTheta + (1.0 / Pe) * laplacianTheta) * dt;
            }
        }

        //----------------------------------------------------------------------------
        // Vorticity Transport Equation
        //----------------------------------------------------------------------------
        OmegaP = Omega;
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                real ddxTheta         = fdOperator::deltaBar_x(ThetaP, i, j)    / (2 * h);
                real ddyPsi           = fdOperator::deltaBar_y(Psi, i, j)       / (2 * h);
                real ddxOmegaP        = fdOperator::deltaBar_x(OmegaP, i, j)    / (2 * h);
                real ddxPsi           = fdOperator::deltaBar_x(Psi, i, j)       / (2 * h);
                real ddyOmegaP        = fdOperator::deltaBar_y(OmegaP, i, j)    / (2 * h);
                real laplacianOmegaP  = (OmegaP[i + 1][j] + OmegaP[i - 1][j] + OmegaP[i][j + 1] + OmegaP[i][j - 1] - 4 * OmegaP[i][j]) / (h * h);
                real buoyancyTerm     = Ri * ddxTheta; 

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
            real error = 0.0;
            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    error = std::max(error, std::abs(Omega[i][j] - OmegaP[i][j]));
                }
            }

                // Print error every 100 iterations
            if (iter % 100 == 0)
                Info << "Iteration: " << iter << " | Error: " << error << std::endl;
            if (error < maxe) {
                Info << "Converged after " << iter << " iterations." << std::endl;
                break;
            }
        }
    }
    //----------------------------------------------------------------------------
    // Create velocities from stream-function obtained above
    //----------------------------------------------------------------------------
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            u[i][j] = fdOperator::deltaBar_y(Psi, i, j) / (2 * h);//(Psi[i][j + 1] - Psi[i][j - 1]) / (2 * h);
            v[i][j] = (-Psi[i + 1][j] + Psi[i - 1][j]) / (2 * h);
        }
    }
    for (int i = 1; i < Nx - 1; ++i) {
        u[i][Ny - 1] = UWall; 
    }

    //----------------------------------------------------------------------------
    // Call VTK Writer
    //----------------------------------------------------------------------------
    writeToVTK(x, y, u, v, Psi, Omega, Theta, "Gr1000Re100cavityResults.vtk");
    return 0;
}