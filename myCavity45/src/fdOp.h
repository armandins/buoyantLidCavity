#ifndef FD_OPERATOR_H
#define FD_OPERATOR_H

#include <iomanip>
#include <vector>

using out = std::ostream&;
out Info = std::cout;
out cErr = std::cerr;

using int16 = int16_t;
using int32 = int32_t;
using int64 = int64_t;
using real  = double;

// Vector 
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
using twoDVec = std::vector<std::vector<double>>;
using oneDVec = std::vector<double>;
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// Modern Array
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <size_t Nx>
using oneDArray = std::array<double, Nx>;

template <size_t Nx, size_t Ny>
using twoDArray = std::array<std::array<double, Ny>, Nx>;
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++



namespace fdOperator {
    // First-order forward difference in x-direction
    double Delta_x(twoDVec& u, int i, int j) {
        return (u[i + 1][j] - u[i][j]);
    }
    // First-order forward difference in y-direction
    double Delta_y(twoDVec& u, int i, int j) {
        return (u[i][j + 1] - u[i][j]);
    }
    
    // First-order backward difference in x-direction
    double Nabla_x(twoDVec& u, int i, int j) {
        return (u[i][j] - u[i - 1][j]);
    }

    // First-order backward difference in y-direction
    double Nabla_y(twoDVec& u, int i, int j) {
        return (u[i][j] - u[i][j - 1]);
    }

    // First-order central difference in x-direction
    double deltaBar_x(twoDVec& u, int i, int j) {
        return (u[i + 1][j] - u[i - 1][j]);
    }

    // First-order central difference in y-direction
    double deltaBar_y(twoDVec& u, int i, int j) {
        return (u[i][j + 1] - u[i][j - 1]);
    }

    // First-order central difference in x-direction
    double delta_x(twoDVec& u, int i, int j) {
        return (u[i + 0.5][j] - u[i - 0.5][j]);
    }

    // First-order central difference in y-direction
    double delta_y(twoDVec& u, int i, int j) {
        return (u[i][j + 0.5] - u[i][j - 0.5]);
    }

    // Second-order central difference in x-direction
    double delta2_x(twoDVec& u, int i, int j) {
        return (u[i + 1][j] - 2 * u[i][j] + u[i-1][j]);
    }

    double delta2_y(twoDVec& u, int i, int j) {
        return (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]);
    }

    // Averaging operator in x-direction
    double mu_x(twoDVec& u, int i, int j) {
        return ((u[i + 0.5][j] + u[i - 0.5][j]) / 2);
    }
    // Averaging operator in y-direction
    double mu_y(twoDVec& u, int i, int j) {
        return ((u[i][j + 0.5] + u[i][j - 0.5]) / 2);
    }
} 

#endif // FD_OPERATOR_H