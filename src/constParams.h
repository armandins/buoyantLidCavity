#ifndef CONST_PARAMS_H
#define CONST_PARAMS_H
#include "fdOp.h"

int16     Nx              { 129           };              // Number of grid points in x-dir
int16     Ny              { 129           };              // Number of grid points in y-dir
int32     maxIt           { 600000        };              // Maximum iter. allowed
real      L               { 1.0           };              // Length of each wall
real      UWall           { 1.0           };              // Upper-wall velocity
real      dt              { 5E-05         };              // Timestep size
real      maxe            { 1E-7          };              // Max-error
real      Re              { 100           };              // Reynolds
real      Pr              { 0.70          };              // Prandtl
real      Gr              { 1E+3          };              // Grashof
real      Pe              { Re*Pr         };              // Peclet
real      Ri              { Gr / (Re * Re)};              // Richardson

#endif

