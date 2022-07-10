#include "ALG.h"

V2D Inv22(const V2D& A)
{
  V2D C; C.resize(2, V1D (2,0.0));
  double Determinante;
  Determinante=(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
  C[0][0]=(+A[1][1])/Determinante;
  C[0][1]=(-A[0][1])/Determinante;
  C[1][0]=(-A[1][0])/Determinante;
  C[1][1]=(+A[0][0])/Determinante;
  return C; 
};

V2D MM(const V2D&  A, const V2D& B)
{
  V2D C;
  C.resize(A.size(), V1D (B[0].size(),0.0));
  for (int i=0; i<A.size(); i++) 
  {
    for (int j=0; j<B[0].size(); j++) 
    {
      for (int k=0; k<A[0].size(); k++) 
      {
        C[i][j]=C[i][j]+A[i][k]*B[k][j];
      }
    }
  }
  return C;
};

V1D MV(const V2D&  M, const V1D& V)
{
  V1D C;
  C.resize(M.size(),0.0);
  for (int i=0; i<M.size(); i++)
  {
    for (int j=0; j<M[0].size(); j++)
    {
      C[i]=C[i]+M[i][j]*V[j];
    }
  }
  return C;
};

V2D MtM(const V2D&  A, const V2D& B)
{
  V2D C;
  C.resize(A[0].size(), V1D (B[0].size(),0.0));
  for (int i=0; i<A[0].size(); i++) 
  {
    for (int j=0; j<B[0].size(); j++) 
    {
      for (int k=0; k<A.size(); k++) 
      {
        C[i][j]=C[i][j]+A[k][i]*B[k][j];
      }
    }
  }
  return C;
};

V2D Mt(const V2D&  A)
{
  V2D C;
  C.resize(A[0].size(), V1D (A.size(),0.0));
  for (int i=0; i<A[0].size(); i++) 
  {
    for (int j=0; j<A.size(); j++) 
    {
      C[i][j]=A[j][i];
    }
  }
  return C;
};

V2D MtpM(const V2D&  A, const V2D& B)
{
  V2D C;
  C.resize(A[0].size(), V1D (A.size(),0.0));
  for (int i=0; i<A[0].size(); i++) 
  {
    for (int j=0; j<A.size(); j++) 
    {
      C[i][j]=A[j][i]+B[i][j];
    }
  }
  return C;
};

double Det( const V2D& m ) 
{
  double res=0.0;
  V1I Ind;
  Ind.resize(2*m.size(),0.0);
  for (int i=0; i<m.size(); i++) 
  {
    Ind[i]=i;
    Ind[m.size()+i]=i;
  }

  int k;
  if (m.size()>2) {k=0;} else {k=1;}
    for (int i=0; i<m.size()-k; i++) 
    {
      double a=1.0;
      double b=1.0;
      for (int j=0; j<m[0].size(); j++) 
      {
        a=a*m[Ind[j]][Ind[j+i]];
        b=b*m[Ind[m.size()-1-j]][Ind[j+i]];
      }
    res=res+a-b;
    }
  return res;
};

double determinant(const V2D vect)
{
  if(vect.size() != vect[0].size())
  {
	throw std::runtime_error("Matrix is not quadratic");
  }
  int dimension = vect.size();

  if(dimension == 0)
  {
	return 1;
  }

  if(dimension == 1)
  {
	return vect[0][0];
  }

  if(dimension == 2)
  {
	return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
  }

  double result = 0;
  int sign = 1;
  for(int i = 0; i < dimension; i++)
  {
	  V2D subVect(dimension - 1, V1D (dimension - 1));
	  for(int m = 1; m < dimension; m++)
	  {
	    int z = 0;
	    for(int n = 0; n < dimension; n++)
	    {
		    if(n != i)
		    {
		    subVect[m-1][z] = vect[m][n];
		    z++;
		    }
	    }
	  }
	  //recursive call
	  result = result + sign * vect[0][i] * determinant(subVect);
	  sign = -sign;
  }
  return result;
}

V2D transpose(const V2D matrix1)
{
  //Transpose-matrix: height = width(matrix), width = height(matrix)
  V2D solution(matrix1[0].size(), V1D (matrix1.size()));
  //Filling solution-matrix
  for(std::size_t i = 0; i < matrix1.size(); i++)
  {
	  for(std::size_t j = 0; j < matrix1[0].size(); j++)
	  {
	    solution[j][i] = matrix1[i][j];
	  }
  }
  return solution;
}

V2D cofactor(const V2D vect)
{
  if(vect.size() != vect[0].size())
  {
	throw std::runtime_error("Matrix is not quadratic");
  }
  V2D solution(vect.size(), V1D (vect.size()));
  V2D subVect(vect.size() - 1, V1D (vect.size() - 1));
  for(std::size_t i = 0; i < vect.size(); i++)
  {
	  for(std::size_t j = 0; j < vect[0].size(); j++)
	  {
	  int p = 0;
	    for(std::size_t x = 0; x < vect.size(); x++)
	    {
		    if(x == i)
		    {
		    continue;
	  	  }
		    int q = 0;
		      for(std::size_t y = 0; y < vect.size(); y++)
		      {
		        if(y == j)
		        {
			      continue;
		        }
		        subVect[p][q] = vect[x][y];
		        q++;
		      }
		    p++;
	    }
	  solution[i][j] = pow(-1, i + j) * determinant(subVect);
	  }
  }
  return solution;
}

V2D inverse(const V2D A)
{
  if(determinant(A) == 0)
  {
	throw std::runtime_error("Determinant is 0");
  }

  double d = 1.0/determinant(A);
  V2D solution(A.size(), V1D (A.size()));
  if(A.size() == 1)
  {
	  V1D ans = {0};
  	ans[0] = 1.0/determinant(A);
  	solution[0] = ans;
	  return solution;
  }
  for(std::size_t i = 0; i < A.size(); i++)
  {
	  for(std::size_t j = 0; j < A.size(); j++)
	  {
	    solution[i][j] = A[i][j];
	  }
  }
  solution = transpose(cofactor(solution));
  for(std::size_t i = 0; i < A.size(); i++)
  {
	  for(std::size_t j = 0; j < A.size(); j++)
	  {
	    solution[i][j] *= d;
	  } 
  }
  return solution;
}

void printMatrix(const V2D vect)
{
  for(std::size_t i = 0; i < vect.size(); i++)
  {
	  for(std::size_t j = 0; j < vect[0].size(); j++)
	  {
	    std::cout << std::setw(8) << vect[i][j] << " ";
	  }
	  std::cout << "\n";
  }
}
