#ifndef HALGH
#define HALGH
#include <vector>   //std::vector
#include <fstream>  //entrada e saida em arquivos
#include <iostream> //std::cout
#include <math.h> // fabs cos sin acos ...
#include <iomanip>
#include <stdexcept>
#include <stdio.h>

typedef std::vector< int > V1I;
typedef std::vector< double > V1D;
typedef std::vector< V1D  > V2D;
typedef std::vector< V2D  > V3D;
typedef std::vector< V3D  > V4D;
V2D Inv22(const V2D&  A);
V2D Inv33(const V2D& A);
V2D MM   (const V2D&  A, const V2D& B);
V1D MV   (const V2D&  M, const V1D& V);
V2D MtM  (const V2D&  A, const V2D& B);
V2D Mt   (const V2D&  A);
V2D MtpM (const V2D&  A, const V2D& B);
double determinant(const V2D vect);
V2D transpose(const V2D matrix1);
V2D cofactor(const V2D vect);
V2D inverse(const V2D vect);
void printMatrix(const V2D vect);
 
class vetor
{
  private:
    int sz;	// número de elementos 
    int* v;	// ponteiro para inteiros 
  public:
    vetor(int);	// construtor 
    ~vetor();	// destrutor
};

double Det( const V2D& m );
#endif
