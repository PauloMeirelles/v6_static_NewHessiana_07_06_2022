#ifndef HQUAH
#define HQUAH

#include <vector>
#include <fstream>
#include <math.h>
#include "ALG.h"

class CQuadratura  //Dunavant/Hammer Triangulo (2D)
{
  private:
	int          Quad;                         ///< Numero da quadratura;
	int    NumberQuad;                         ///< Numero de pontos da quadratura;
	int  DegreeApprox;                         ///< Grau de aproximacao da quadratura;
	int          nXsi;                         ///< Numero de coordenadas adimensionais Qsi;
	V2D           Xsi;                         ///< Vetor 2D de Coordenadas Qsi;
	V1D             W;                         ///< Pesos da quadratura;
	void FPDunavant();                         ///< Funcao que gera os pontos de Dunavant/Hammer no intervalo de 0 a 1 no espa�o homogeneo para (pontos de integracao 2D) grau 1 a 20;
  public:
	CQuadratura(int Quad, int& DegreeApprox);        ///< Funcao construtora das coordenadas Xsi e dos respectivos pesos da Quadratura;
	~CQuadratura();                                  ///< Destrutor da classe Quadratura
	int r_Quad () {return (Quad);}                   ///< Numero da quadratura
	int r_NumberQuad () {return (NumberQuad);}       ///< Numero de pontos da quadratura;
	int r_DegreeApprox () {return (DegreeApprox);}   ///< Grau de aproximacao da quadratura;
	int r_nXsi () {return (nXsi);}                   ///< Numero de coordenadas adimensionais xsi;
	V1D r_Xsi (int pto) {return (Xsi[pto]);}         ///< Retorna a coordenada de Xsi(i);
	double r_W (int i) {return (W[i]);}              ///< Retorna o peso de gauss W(i);
};

typedef std::vector< CQuadratura> tvQuadratura; 	//Vetor de Quadratura;

void Le_Quadraturas(tvQuadratura& Qua, const std::string& NAr);   //Verificar se NAr uma var para leitura de entrada de dados

#endif
