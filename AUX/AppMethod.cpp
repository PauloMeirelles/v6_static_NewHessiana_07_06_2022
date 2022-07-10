#include "AppMethod.h"

APPMETHOD::APPMETHOD(const int& _Method, const int& _OutputFile)
{
  Method=_Method;
  OutputFile=_OutputFile;
}

void Le_AppMethod(tvAPPMETHOD &cAppMethod, const std::string& NAr)
{
  int nMethod;         ///< Numero do metodo de aproximação pos process
  int nOutputFile;     ///< Numero do modelo do arquivo de saida
  char s[1000];        ///< Variavel auxiliar
  std::ifstream Ent;   ///< Arquivo de leirura
  Ent.open(NAr.c_str());
  Ent.getline(s,1000);
  Ent >>  nMethod ; Ent.getline(s,1000);
  Ent.getline(s,1000);
  Ent >>  nOutputFile ; Ent.getline(s,1000);
  APPMETHOD AppMethod1(nMethod, nOutputFile);
  cAppMethod.push_back(AppMethod1);  
  Ent.close();
}

APPMETHOD::~APPMETHOD()
{
}
