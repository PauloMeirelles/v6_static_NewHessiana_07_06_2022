#include "MAT.h"

CMa::CMa(const int& _NMat, const double& _E0, const double& _Nu0, const double& _rho, const double& _Mi)
{
  NMat=_NMat;
  E0=_E0;
  Nu0=_Nu0;
  Et=E0/(2.0*(1.0+Nu0));
  Kd=E0/(3.0*(1.0-2.0*Nu0));
  LameM=Et;
  LameL=2.0*Et*Nu0/(1.0-2.0*Nu0);
  rho=_rho;
  Mi=_Mi;
}

tvMa Le_Mat(const std::string& NAr)
{
  tvMa Ma;
  int i,nMa;
  char s[1000];
  std::ifstream Ent; // para leitura
  Ent.open(NAr.c_str());
  Ent.getline(s,1000);
  Ent >>  nMa; Ent.getline(s,1000);
  Ent.getline(s,1000);

  for (i=0; i<nMa; i++)
  {
	  double _E0;
	  double _Nu0;
	  double _rho;
	  double _Mi;
	  int NMat;
	  Ent >> NMat >> _E0 >> _Nu0 >> _rho >> _Mi; Ent.getline(s,1000);
	  Ma.push_back(CMa(NMat,_E0,_Nu0,_rho,_Mi));
  }

  Ent.close();
  return Ma;
}
