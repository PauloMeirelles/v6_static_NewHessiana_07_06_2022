#ifndef HAPPMETHODH
#define HAPPMETHODH

#include "ALG.h"
#include "NOS.h"
#include "MAT.h"
#include "QUA.h"
#include "CCO.h"
#include "solution.h"
#include <vector>   //std::vector
#include <fstream>  //entrada e saida em arquivos
#include <iostream> //std::cout

class APPMETHOD
{
  private:
    int Method;              //> Metodo de aproximação pós-processamento
    int OutputFile;          //> Arquivo de saida em qual formato
  public: 
    int r_Method(){return Method;}
    int r_OutputFile(){return OutputFile;}

    APPMETHOD(const int& _Method, const int& _OutputFile);
    ~APPMETHOD();
};

typedef std::vector< APPMETHOD > tvAPPMETHOD; //Vetor de Elementos finitos

    
void Le_AppMethod(tvAPPMETHOD &cAppMethod, const std::string& NAr);

#endif