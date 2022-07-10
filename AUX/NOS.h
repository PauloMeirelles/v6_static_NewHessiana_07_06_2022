#ifndef HNOSH       // Equivalente a "#if !defined" -> verificam somente a presença ou ausência de identificadores definidos com #define
#define HNOSH       // Idenfiticador do cabeçalho NOS

#include <vector>   ///<std:: Para vetorização;
#include <fstream>  ///<std:: Para leitura de arquivos;
#include <iostream> ///<std:: Cout and cin;
#include "ALG.h"    ///<Inclui o cabeçalho ALG.h;

///Classe: Nos Geometricos Elemento Finito Triangular (2D);
class CNo
{
  private:
	int            NumberNo;   ///< Numero do No;
	double          X[2][2];   ///< Coordenadas iniciais(Configuracao 0 [X_{1} , X_{2}]) e finais(Configuração 1 [Y_{1} , Y_{2}]) do No;
	int     CondCountour[2];   ///< Condicao de contorno (1->vinculado, 0->livre) ;
	double vCondCountour[2];   ///< Condicao de contorno (valor do deslocamento prescrito)   ;
	double PConcentrated[2];   ///< Carga concentrada no nó;
	int           Adress[2];   ///< Enderecamento das coordenadas globais (Utilizado para definir na matriz global os nó que farão parte ou não);
	int                  FP;   ///< 0:Não faz parte da estrutura ou 1:Faz parte da estrutura;
	int 		  dimension;   ///< DImensão do problema

	//-------------------------/// Verificar se FP deve fazer parte do código nesse momento;
	V1D    no_E;               ///< Deformação de Green-Lagrange do cada nó (vetor de 4 termos);
	V1D    no_S;               ///< Tensão de PiolaKirchhoff de Segunda Espécie de cada nó (vetor de 4 termos -> Caso 2D);
	V1D    no_Sig;             ///< Tensão de Cauchy
	double no_ue;              ///< Energia específica de deformação de cada nó ;
	V1I    Elms;               ///< Armazena uma lista com os elementos que contém este nó;
	V1I    no_Elms;            ///< Armazena uma lista com número local no no nos elementos que contém este nó;

  public:
    int    r_dimension() {return(dimension);} 
	int    r_NumberNo    () {return(NumberNo);}                                       ///<Retorna (Le) o numero do No;
	double r_X    (int estado, int dir) {return(X[estado][dir]);}                     ///<Retorna (Le) a coordenada X(i) do No no estado 0 ou 1;
	void   w_X    (int estado, int dir, double nX) {X[estado][dir]=nX;}               ///<Altera (Escreve) a coordenada X(i) do No no estado 0 ou 1;
	int    r_CondCountour   (int dir)             {return (CondCountour[dir]);}       ///<Retorna (Le) condicoes de contorno;
 	void   w_CondCountour   (int dir, int value)    {CondCountour[dir]=value;}        ///<Altera (Escreve) condicoes de contorno;
	double r_vCondCountour  (int dir)             {return (vCondCountour[dir]);}      ///<Retorna (Le) valor das condicoes de contorno prescritas;
	void   w_vCondCountour  (int dir, double value) {vCondCountour[dir]=value;}       ///<Altera (Escreve) valor das condicoes de contorno prescritas;
	int    r_Adress   (int dir) {return (Adress[dir]);}                               ///<Retorna (Le) endereçamento dos Nos em relacao a matriz global;
	void   w_Adress   (int dir, int value) {Adress[dir]=value;}                       ///<Altera (Escreve) endereçamento dos Nos em relacao a matriz global;
	double r_PConcentrated    (int dir) {return (PConcentrated[dir]);}                ///<Retorna (Le) valor das cargas concentradas nos Nos;
	void   w_PConcentrated    (int dir, double value) {PConcentrated[dir]=value;}     ///<Altera (Escreve) valor das cargas concentradas nos Nos;
    // int    r_StepForce	() {return (StepForce);}
    // void   w_StepForce	(int value){StepForce=value;}	
	double r_FP   (       ) {return (FP);}                                            ///<Retorna (Le) se o No faz parte ou nao da estrutura;
	void   w_FP   (int value) {FP=value;}                                             ///<Altera (Escreve) o No que faz parte ou nao da estrutura;
	double r_no_E (int k1 ) {return (no_E[k1]);}                                      ///<Retorna (Le) as Deformacoes de Green do No;
	void   w_no_E (int k1, double Value) {no_E[k1]=Value;}                            ///<Altera (Escreve) as Deformacoes de Green do No;
	double r_no_S (int k1 ) {return (no_S[k1]);}                                      ///<Retorna (Le) as Tensoes de PK no No;
	void   w_no_S (int k1, double Value) {no_S[k1]=Value;}                            ///<Altera (Escreve) as Tensoes de PK no No;
	double r_no_Sig (int k1 ) {return (no_Sig[k1]);}                                      
	void   w_no_Sig (int k1, double Value) {no_Sig[k1]=Value;}                            
	double r_no_ue() {return (no_ue);}                                                ///<Retorna (Le) a energia especifica de deformacao do No;
	void   w_no_ue(double Value) {no_ue=Value;}                                       ///<Altera (Escreve) a energia especifica de deformacao do No;
	void   I_Elm(int Elm, int NoLocal) { Elms.push_back(Elm); no_Elms.push_back(NoLocal); }
	int    r_Elm(int k1) {return Elms[k1]; }
	int    r_no_Elm(int k1) {return no_Elms[k1]; }
    int    sizeElm() {return Elms.size(); }
	CNo(const int& , const double& , const double& , const double& , const double&);  ///<Funcao contrtutora;
	~CNo();                    ///< Classe Destrutora CNo ;
};

typedef std::vector< CNo > tvNo;	//Vetor Nos Geometricos
void Le_Nos(tvNo& No, const std::string& NAr);
void Le_CC(tvNo& No, const std::string& NAr, int& ngl);
void Le_P(tvNo& No, const std::string& NAr);


#endif
