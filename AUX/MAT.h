#ifndef HMATH
#define HMATH

#include <vector>
#include <fstream>
#include "ALG.h"

///Classe: Material
class CMa
{
  private:
	int     NMat;       ///<Numero do Material
	double    E0;	    ///<Modulo de elasticidade
	double   Nu0;	    ///<Coeficiente de Poisson
	double    Et;	    ///<Modulo de elasticidade trasnversal
	double    Kd;	    ///<Bulk Modulus (Módulo volumétrico)
	double LameM;	    ///<Segunda Constante de lame = μ = G = Et
	double LameL;	    ///<Primeira Constante de Lame = λ = 2*G*Nu/(1-2*Nu)=E*Nu/((1+Nu)*(1-2*Nu))
	double   rho;       ///<Massa especifica
	double    Mi;       ///<Viscodidade dinâmica
  public:
	int    r_NMat() {return (NMat);}      ///<Escreve a numera��o do Material;
	double r_E0() {return (E0);}          ///<Escreve o E0;
	double r_Nu0() {return (Nu0);}        ///<Escreve o Nu0;
	double r_Et() {return (Et);}          ///<Escreve o Et;
	double r_Kd() {return (Kd);}          ///<Escreve o Kd;
	double r_LameM() {return (LameM);}    ///<Escreve o LameM;
	double r_LameL() {return (LameL);}    ///<Escreve o LameL;
	double r_rho() {return (rho);}        ///<Escreve o Rho;
	double r_Mi() {return (Mi);}          ///<Escreve o Mi;
	CMa(const int& NMat, const double& E0, const double& Nu0, const double& rho, const double& Mi);
};
typedef std::vector< CMa > tvMa;	//Vetor Material
tvMa  Le_Mat(const std::string& NAr);

#endif