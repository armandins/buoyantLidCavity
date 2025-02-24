#ifndef VTKWRITER_H
#define VTKWRITER_H

#include <vector>
#include <string>
#include <fstream>
#include <fstream>
#include "fdOp.h"

// ===============================================================================
// VTK Format Writer (*.vtk)                                                       
// ===============================================================================

void writeToVTK(const oneDVec& x, const oneDVec& y, 
            const twoDVec& u, 
            const twoDVec& v, 
            const twoDVec& Psi, 
            const twoDVec& Omega,
            const twoDVec& Theta,  
            const std::string& filename) {
std::ofstream vtkFile(filename);
if (!vtkFile.is_open()) {
    cErr << "Error: Could not open file " << filename << " for writing." << '\n';
    return;
}

int Nx = x.size();
int Ny = y.size();

vtkFile << "# vtk DataFile Version 3.0\n";
vtkFile << "Lid-Driven Cavity Flow\n";
vtkFile << "ASCII\n";
vtkFile << "DATASET STRUCTURED_GRID\n";
vtkFile << "DIMENSIONS " << Nx << " " << Ny << " 1\n";
vtkFile << "POINTS " << Nx * Ny << " double\n";

for (int j = 0; j < Ny; ++j) {
    for (int i = 0; i < Nx; ++i) {
        vtkFile << x[i] << " " << y[j] << " 0.0\n";
    }
}

vtkFile << "POINT_DATA " << Nx * Ny << "\n";
vtkFile << "VECTORS velocity double\n";
for (int j = 0; j < Ny; ++j) {
    for (int i = 0; i < Nx; ++i) {
        vtkFile << u[i][j] << " " << v[i][j] << " 0.0\n";
    }
}

vtkFile << "SCALARS stream_function double 1\n";
vtkFile << "LOOKUP_TABLE default\n";
for (int j = 0; j < Ny; ++j) {
    for (int i = 0; i < Nx; ++i) {
        vtkFile << Psi[i][j] << '\n';
    }
}

vtkFile << "SCALARS u_velocity double 1\n";
vtkFile << "LOOKUP_TABLE default\n";
for (int j = 0; j < Ny; ++j) {
    for (int i = 0; i < Nx; ++i) {
        vtkFile << u[i][j] << '\n';
    }
}

vtkFile << "SCALARS v_velocity double 1\n";
vtkFile << "LOOKUP_TABLE default\n";
for (int j = 0; j < Ny; ++j) {
    for (int i = 0; i < Nx; ++i) {
        vtkFile << v[i][j] << '\n';
    }
}

vtkFile << "SCALARS vorticity double 1\n";
vtkFile << "LOOKUP_TABLE default\n";
for (int j = 0; j < Ny; ++j) {
    for (int i = 0; i < Nx; ++i) {
        vtkFile << Omega[i][j] << '\n';
    }
}

vtkFile << "SCALARS theta double 1\n";
vtkFile << "LOOKUP_TABLE default\n";
for (int j = 0; j < Ny; ++j) {
    for (int i = 0; i < Nx; ++i) {
        vtkFile << Theta[i][j] << '\n';
    }
}

vtkFile.close();
Info << "Results written to " << filename << '\n';
}

#endif // VTKWRITER_H
