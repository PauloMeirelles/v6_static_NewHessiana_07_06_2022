#ifndef HTEMPOH
#define HTEMPOH

#include "ALG.h"
#include "NOS.h"
#include "MAT.h"
#include "QUA.h"
#include "CCO.h"
#include "solution.h"
#include <vector>   //std::vector
#include <fstream>  //entrada e saida em arquivos
#include <iostream> //std::cout

class TEMPO
{
  private:
    int TimeStep;            //> Passo de Tempo
    int NTimeStep;           //> Numero de Passos de Tempo
    int ForceStep;           //> Numero de incrementos de passo de força
    double Beta, Gamma;      //> Gamma e Beta para Newmark
    double Tolerance;        //> Tolerancia NR com controle de posição

  public: 
    int r_TimeStep(){return TimeStep;}
    int r_NTimeStep(){return NTimeStep;}
    int r_ForceStep(){return ForceStep;}
    double r_Tolerance(){return Tolerance;}
    double r_Beta(){return Beta;}
    double r_Gamma(){return Gamma;}
    TEMPO(const int& _TimeStep, const int& _NTimeStep,  const int& _ForceStep, const double& _Tolerance, const double& _Beta,  const double& _Gamma);
    ~TEMPO();
};

typedef std::vector< TEMPO > tvTempo; //Vetor de Elementos finitos

    
void Le_Tempo(tvTempo &vTempo, const std::string& NAr);

#endif
