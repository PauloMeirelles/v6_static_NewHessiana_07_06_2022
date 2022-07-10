#ifndef HMUMPSSOLVEH
#define HMUMPSSOLVEH

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <petsc.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscviewer.h>
#include <vector>
#include "Eigen/Dense"
#include <Eigen/Core> 
#include <Eigen/SVD>
#include "ALG.h"

V1D MUMPSsolver(int argc, char *argv[], V2D& Ma, V1D& Vb);

#endif