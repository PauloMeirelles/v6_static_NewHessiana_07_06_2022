#include "TEMPO.h"

TEMPO::TEMPO(const int& _TimeStep, const int& _NTimeStep,  const int& _ForceStep, const double& _Tolerance, const double& _Beta,  const double& _Gamma)
{
  TimeStep = _TimeStep;
  NTimeStep = _NTimeStep;
  Beta = _Beta;
  Gamma = _Gamma;
  ForceStep = _ForceStep;
  Tolerance = _Tolerance;
}

void Le_Tempo(tvTempo &vTempo, const std::string& NAr)
{
  int nTimeStep=0;       ///< Numero de passos de tempo
  int TimeStep=0;        ///< Valor do passos de tempo
  int ForceStep=0;       ///< Valor do passos de tempo
  double Beta, Gamma, Tolerance; Beta=0.0; Gamma=0.0; Tolerance=0.0;  ///< Beta e Gamma Newmark Beta, Tolerancia
  char s[1000];        ///< Variavel auxiliar
  std::ifstream Ent;   ///< Arquivo de leirura
  Ent.open(NAr.c_str());
  Ent.getline(s,1000);
  Ent >>  ForceStep ; Ent.getline(s,1000);
  Ent.getline(s,1000);
  Ent >>  Tolerance ; Ent.getline(s,1000);
  Ent.getline(s,1000);
  Ent >>  nTimeStep ; Ent.getline(s,1000);
  Ent.getline(s,1000);
  Ent >> TimeStep >> Beta >> Gamma; Ent.getline(s,1000);
  TEMPO tempo(TimeStep,nTimeStep,ForceStep,Tolerance,Beta,Gamma);
  vTempo.push_back(tempo);  
  Ent.close();
}

TEMPO::~TEMPO()
{
}


