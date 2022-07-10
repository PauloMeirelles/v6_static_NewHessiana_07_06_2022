#ifndef HELMH
#define HELMH

#include "ALG.h"
#include "NOS.h"
#include "MAT.h"
#include "QUA.h"
#include "CCO.h"
#include "solution.h"
#include <vector>   //std::vector
#include <fstream>  //entrada e saida em arquivos
#include <iostream> //std::cout

class CEl
{
  private:
	int NEl;                    ///<Numero do elemento finito
	int NNo;                    ///<Numero de nos
	int NNEl1m;                 ///<numero de nos do elemento com grau do polinomio uma vez menor
	int dimension;              ///<Dimensoes do elemento 1:1D 2:2D ou 3:3D;
	int DegreeApproxPol;        ///<Grau do polinomio aproximador do elemento;
	int dDegreeApproxPol;       ///<Grau da primeira derivada do polinomio aproximador do elemento;
	std::vector< CNo * > No;	///<Nos do elemento(ponteiro)
	std::vector< int > no;      ///<Numero global dos nos do elemento
	CMa* Ma;                    ///<Acesso a classe do Material do elemento (ponteiro)
	int NMa;                    ///<Numero do Material do elemento
	CQuadratura *Quad;          ///<Acesso a classe de quadratura de integracao (ponteiro)
	int NQuad;                  ///<Numero da Quadratura de integracao
	CCo* Config;                ///<Acesso configuracoes como arquivo de entrada (ponteiro)

	// Variáveis armazenadas para cada nó do elemento:
	std::vector< V1D > Xsi;     ///<Coordenadas adimensionais dos nós do elemento no espaço homogêneo
	V1D Fint;                   ///< Vetor de forças internas
	V2D Hessiana;               ///< Matrix Hessiana
	std::vector< V1D >  no_E;   ///< Tensor de Deformação de Green-Lagrange
	std::vector< V1D >  no_S;   ///< Tensor de Tensões de PiolaKirchhoff de segunda espécie
	std::vector< V1D >  no_Sig; ///< Tensor de Tensões de PiolaKirchhoff de segunda espécie
	std::vector< double > no_ue;///< Energia específica de deformação
	V2D Mass;                   ///< Matrix de Massa
	V2D MRP;                    ///< MRP

	// Variáveis armazenadas por elemento:
	double Ue;                  ///< Energia de deformação
	// Variáveis armazenadas para cada ponto de integração do elemento:
	std::vector< V1D >  PHI;    ///< Valor das funções de forma nos pontos de integração
	std::vector< V2D >  DPHI;   ///< Valor da primeira derivada das funcoes de forma nos pontos de integração
	std::vector< V1D >  X;      ///< Posição do ponto de integração na configuração inicial
	std::vector< V2D >  A0;     ///< Gradiente da função mudança de configuração (inicial)
	std::vector< double > J0;   ///< Determinante do Gradiente de tranformação inicial (Jacobiano da transformação)
	std::vector< V2D >  invA0;  ///< Matriz inversa da Gradiente de tranformação inicial
    std::vector< V1D >  Y;      ///< Posição do ponto de integração na configuração final
	std::vector< V2D >  A1;     ///< Gradiente de tranformação final
    std::vector< double > J1;   ///< Determinante do Gradiente de tranformação atual (jacobiano da transformação)
	std::vector< V2D >  A;      ///< Gradiente de transformação total (inicial-linha 1 ; final-linha 2)
	std::vector< double > J;    ///< Determinante do Gradiente de tranformação (jacobiano da transformação)
	std::vector< V1D > Sig;     ///< Tensor de Tensões de CAUCHY;
	std::vector< V2D >  C;      ///< Tensor de alongamento ou estiramento de Cauchy-Green
	std::vector< V1D >  E;      ///< Tensor de Deformação de Green-Lagrange;
	std::vector< V2D >  CC;     ///< Tensor constitutivo elástico de quarta ordem
	std::vector< V1D >  S;      ///< Tensor de Tensão de Piola Kirchhoff  de segunda espécie
	std::vector< double > ue;   ///< Energia específica de deformação
	std::vector< V1D >  fint;   ///< Vetor de forças específica* internas
	std::vector< V2D >  h;      ///< Matrix Hesiana específica*
	std::vector< V2D > m;       ///< Matrix de Massa espefífica
	std::vector< V2D > mrp;     ///< Matrix de interpolacao via Residuos Ponderados de variáveis
	std::vector< V1D > vrp;     ///< Vetor de interpolacao via Residuos Ponderados de variáveis

	int c_NNEl1m();

  public:
	int  r_NEl() {return (NEl);}                             ///<Retorna o Numero do elemento finito
	int  r_NNo() {return (NNo);}                             ///<Retorna o Numero de nos
	int  r_NNEl1m() {return (NNEl1m);}                       ///<Retorna o numero de nos do elemento com grau do polinomio uma vez menor
	int  r_dimension() {return (dimension);}                 ///<Retorna a Dimensão do elemento 1 2 ou 3;
	int  r_DegreeApproxPol() {return (DegreeApproxPol);}     ///<Retorna o Grau do polinomio aproximador do elemento;
	int  r_dDegreeApproxPol() {return (dDegreeApproxPol);}   ///<Retorna o Grau da primeira derivada do polinomio aproximador do elemento;
	std::vector< V1D > r_Xsi() {return (Xsi);}               ///<Retorna coordenadas homogeneas na forma alocados em um vetor;

	CNo *r_No(int i) {return (No[i]);}      ///<Retorna os nós geometricos locais(ponteiro)
	int r_no(int i) {return (no[i]);}       ///<Retorna o numero global do no geometrico
	CMa *r_Ma() {return (Ma);}              ///<Retorna o material do elemento (ponteiro)
	int r_NMa() {return (NMa);}           ///<Retorna o número do material do elemento
	CQuadratura *r_Quad() {return (Quad);}  ///<Retorna a Quadratura de integracao (ponteiro)
	int r_NQuad() {return (NQuad);}         ///<Retorna o numero da Quadratura de integracao
	CCo *r_Co() {return (Config);}              ///<Retorna Configuracoes como arquivo de entrada (ponteiro)
	double r_Hessiana(int i, int j) {return (Hessiana[i][j]);}
	double r_Mass(int i, int j) {return (Mass[i][j]);}
    double r_E(int n, int i) {return (E[n][i]);}
    double r_S(int n, int i) {return (S[n][i]);}
	double r_Sig(int n, int i) {return (Sig[n][i]);}
    double r_ue(int n) {return (ue[n]);}
    double r_Fint(int i) {return (Fint[i]);}
    double r_PHI(int iQua,int iNo) {return (PHI[iQua][iNo]); }
	V1D  r_no_E(int nno) { return no_E[nno]; }                       ///< Retorna Deformação de Green-Lagrange no No escolhido
	V1D  r_no_S(int nno) { return no_S[nno]; }                       ///< Retorna o Tensor de Tensões de PiolaKirchhoff
	V1D  r_no_Sig(int nno) { return no_Sig[nno]; }                       ///< Retorna o Tensor de Tensões de PiolaKirchhoff
	void  w_no_E(int nno,int i,double val) {no_E[nno][i]=val; }      ///< Retorna alocando a Deformação de Green-Lagrange no No escolhido
	void  w_no_S(int nno,int i,double val) {no_S[nno][i]=val; }      ///< Retorna alocando o Tensor de Tensões de PiolaKirchhoff
	double r_no_ue(int nno) { return no_ue[nno]; }                   ///< Retorna a Energia específica de deformação no no escolhido
    void   w_no_ue(int nno,double val) {no_ue[nno]=val; }
	double r_no_Y(int nno,int dir) { return No[nno]->r_X(1,dir); }   ///< Retorna a posicao deslocada do no escolhido
	double r_no_U(int nno,int dir) { return (No[nno]->r_X(1,dir)-No[nno]->r_X(0,dir)); }   ///< Retorna o deslocamento do no escolhido

	CEl(const int&, const int&, const int&, tvMa&, const int&, const V1I&, tvNo&, const int&, tvQuadratura&, CCo&);
	void ite();

	void Imp_Rel_Geral();
	std::vector < std::vector <int > >  TriL(const int& DegreeApproxPol); // USO ISSO ?
	~CEl(); ///< Destrutor CEl
};

typedef std::vector< CEl > tvEl; //Vetor de Elementos finitos

V2D c_Xsi(const int& NNo, const int& DegreeApproxPol, const int& dimension);   ///<calcula as Coordenadas adimensionais dos nós do elemento no espaço homogêneo
V1D calc_PHI (const int& DegreeApproxPol, const int& NNo, const int dimension, const V1D& Xsi/*, const V2D&  FF*/); ///< Calcula o Valor das funções de forma nos pontos de integração
V2D calc_DPHI(const int& DegreeApproxPol, const int& NNo, const int dimension, const V1D& Xsi/*, const V3D& dFF*/); ///< Calcula o Valor da primeira derivada das funcoes de forma nos pontos de integração
V1D calc_Xn  (const int& dimension, const V1D&  PHI, std::vector< CNo * > No, const int estado);  ///< Calcula a Posição do ponto de integração na configuração inicial
V2D calc_An  (const V2D& DPHI, std::vector< CNo * > No, const int estado);  ///< Calcula o Gradiente de tranformação inicial ou final (n=estado)
V1D calc_Ei  (const V2D& Ci);
V2D calc_CC  (const int& dimension, CMa *Ma);
V1D calc_S   (const int& dimension, const V1D& Ei, V2D& CCi);
V1D calc_Sig (const int& dimension, const V2D& A,const int& J,const V1D& S);
double c_ue(const V1D& E, CMa *Ma);

void calc_Fint_H_M(const V1D& PHI, const V2D& DPHI, const V2D& Di, const V2D& Ai, const V1D& Si, const V2D& CCi, V1D& finti, V2D& hi, V2D& mi, V2D& mrpi, V1D& vrpi, V2D const A1);

void cpto(const int& DegreeApproxPol,const int& NNo,const int& dimension,const V1D XsiNo, std::vector< CNo * > No, CMa *Ma, V1D& E, V1D& S, V1D& Sig, double& ue);

void Le_Elementos(tvEl& El, tvNo& No, tvMa& Ma, tvQuadratura& Qua, CCo& Co);


#endif
