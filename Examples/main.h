//  Created by Wilson Wesley Wutzow and adapted by Paulo Henrique de Freitas Meirelles on 09/12/11.
//  Copyright 2011 Wutzow. All rights reserved.
//
#ifndef HmainH
#define HmainH
//--------------- Bibliotecas ----------------------------
//#include "Eigen/Dense"
#include <time.h>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <string>
#include <cstring>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <errno.h>
#include <stdlib.h>
#include <petsc.h>
#include <petscksp.h>
#include "CCO.h"
#include "ALG.h"
#include "NOS.h"
#include "ELM.h"
#include "MAT.h"
#include "QUA.h"
#include "TEMPO.h"
#include "AppMethod.h"
#include "solution.h"
#include "MumpsSolve.h"

void Gera_VTK_69(CCo& Co, tvNo& No, tvEl& El, const int& ite);    //> Gerador de saida VTK Interpolado
void Gera_VTU_69(CCo& Co, tvNo& No, tvEl& El, const int& ite);    //> Gerador de saida VTK Nao Intepolado
void Gera_VTK_69_NI(CCo &Co, tvNo &No, tvEl &El, const int &ite); //> Gerador de saida VTU Intepolado

#endif
