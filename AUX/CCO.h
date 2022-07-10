#ifndef HCCOH
#define HCCOH

#include <fstream>
#include <iostream> //std::cout
#include <sstream>
#include <string>

///Classe: Configuracao
class CCo
{
  private:
    int  Nconfig;         ///<Numero
    std::string DIR;	    ///<Diretorio
    std::string ArqTEMPO;	///<Nome Arq Tempo - Newmark Beta;
    std::string ArqMETHOD;///<Nome Arq Tempo - Newmark Beta;
    std::string ArqMAT;	  ///<Nome Arq Materiais;
    std::string ArqNOS;	  ///<Nome Arq Nos;
    std::string ArqELM;	  ///<Nome Arq Elementos;
    std::string ArqQUA;	  ///<Nome Arq Quadraturas;
    std::string ArqREL;	  ///<Nome Arq Relatorios;
    std::string ArqVTK;	  ///<Nome Arq VTK;
    std::string ArqVTU;	  ///<Nome Arq VTU;
    std::string ArqCC ;   ///<Nome Arq Condiçoes de Contorno;
    std::string ArqP  ;   ///<Nome Arq Cargas Concentradas;

  public:
    int  r_Nconfig() {return (Nconfig);}        ///<Numero
    std::string r_DIR() {return (DIR);}         ///<Diretorio do exemplo
    std::string r_ArqTEMPO() {return (ArqTEMPO);} ///<Nome Arq Tempo - Newmark Beta;
    std::string r_ArqMETHOD() {return (ArqMETHOD);} ///<Nome Arq Tempo - Newmark Beta;
    std::string r_ArqMAT() {return (ArqMAT);} 	///<Nome Arq Materiais;
    std::string r_ArqNOS() {return (ArqNOS);} 	///<Nome Arq Nos;
    std::string r_ArqELM() {return (ArqELM);} 	///<Nome Arq Elementos;
    std::string r_ArqQUA() {return (ArqQUA);} 	///<Nome Arq Quadraturas;
    std::string r_ArqREL() {return (ArqREL);} 	///<Nome Arq REL;
    std::string r_ArqVTK() {return (ArqVTK);} 	///<Nome Arq VTK;
    std::string r_ArqVTU() {return (ArqVTU);} 	///<Nome Arq VTK;
    std::string r_ArqCC()  {return (ArqCC );} 	///<Nome Arq Condicoes de Contorno;
    std::string r_ArqP()   {return (ArqP  );} 	///<Nome Arq Cargas Concentradas;
    CCo(const std::string& ArqCON );            ///<Construtor da Classe configuração
};

std::string intToString(int i);

#endif
