#ifndef HSOLUTIONH
#define HSOLUTIONH

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
#include "NOS.h"
#include "MAT.h"
#include "QUA.h"
#include "CCO.h"
#include "ELM.h"
#include <fstream>  //entrada e saida em arquivos


V1D KSPresolve(int argc, char *argv[], V2D& Ma, V1D& Vb);
V1D MUMPSresolve(int argc, char *argv[], V2D& Ma, V1D& Vb);
V2D invEIGEN(V2D& Ma);
V2D invAnyEIGEN(V2D& Ma);
V2D pseudoInverse(V2D& Ma); 

V1D PETScSolver(int argc, char *argv[], Mat& H_static,V1D&  g_);
V1D EIGENsolver(int argc, char *argv[], V2D& Hess, V1D& _g);
//V1D StaticEq(int argc, char *argv[], tvEl& _El, const int& ngl, const double& nFS, const double &FStep);

#endif
