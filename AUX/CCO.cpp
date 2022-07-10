#include "CCO.h"

std::string intToString(int i)
{
  std::ostringstream ss;
  std::string s;
  ss << i;
  s = ss.str();
  return s;
}

CCo::CCo(const std::string& ArqCON)
{
  char s[1000];
  std::ifstream ENT; // para leituraF
  ENT.open(ArqCON.c_str());
  ENT.getline(s,1000);
  std::cout << s << std::endl;
  ENT >> DIR; ENT.getline(s,1000);
  ENT.getline(s,1000);
  ENT >> ArqTEMPO; ENT.getline(s,1000);
  ENT.getline(s,1000);
  ENT >> ArqMETHOD; ENT.getline(s,1000);
  ENT.getline(s,1000);
  ENT >> ArqMAT; ENT.getline(s,1000);
  ENT.getline(s,1000);
  ENT >> ArqNOS; ENT.getline(s,1000);
  ENT.getline(s,1000);
  ENT >> ArqCC;  ENT.getline(s,1000);
  ENT.getline(s,1000);
  ENT >> ArqP;   ENT.getline(s,1000);
  ENT.getline(s,1000);
  ENT >> ArqELM; ENT.getline(s,1000);
  ENT.getline(s,1000);
  ENT >> ArqQUA; ENT.getline(s,1000);
  ENT.getline(s,1000);
  ENT >> ArqREL; ENT.getline(s,1000);    
  ENT.getline(s,1000);
  ENT >> ArqVTK; ENT.getline(s,1000);
  ENT.getline(s,1000);
  ENT >> ArqVTU; ENT.getline(s,1000); 
  ArqTEMPO.insert(0,DIR);
  ArqMETHOD.insert(0,DIR);
  ArqMAT.insert(0,DIR);
  ArqNOS.insert(0,DIR);
  ArqCC.insert(0,DIR);
  ArqELM.insert(0,DIR);
  ArqQUA.insert(0,DIR);
  ArqREL.insert(0,DIR);
  ArqVTK.insert(0,DIR);
  ArqP.insert(0,DIR);
  ArqVTU.insert(0,DIR);
  std::cout << "Arquivo de Config.:    ";
  std::cout << ArqCON.c_str() << std::endl;
  std::cout << "Diretorio de Calculo:  ";
  std::cout << DIR.c_str() << std::endl;
  std::cout << "Arquivo de Materiais:   ";
  std::cout << ArqMAT.c_str() << std::endl;
  std::cout << "Arquivo para Sistema Dinâmico (Newmark Beta):   ";
  std::cout << ArqTEMPO.c_str() << std::endl;
  std::cout << "Arquivo para definição dos metodos de PP e Output:   ";
  std::cout << ArqMETHOD.c_str() << std::endl;
  std::cout << "Arquivo de Nos:         ";
  std::cout << ArqNOS.c_str() << std::endl;
  std::cout << "Arquivo de CondCon:     ";
  std::cout << ArqCC.c_str() << std::endl;
  std::cout << "Arquivo de P concentr.: ";
  std::cout << ArqP.c_str()  << std::endl;
  std::cout << "Arquivo de Elementos:   ";
  std::cout << ArqELM.c_str() << std::endl;
  std::cout << "Arquivo de Quadraturas: ";
  std::cout << ArqQUA.c_str() << std::endl;
  std::cout << "Arquivo de Relatorios: ";
  std::cout << ArqREL.c_str() << std::endl;
  std::cout << "Arquivo para o Paraview: ";
  std::cout << ArqVTK.c_str() << std::endl;
  std::cout << "Arquivo Visualization Toolkit Unstructured Points Data (VTU): ";
  std::cout << ArqVTU.c_str() << std::endl;
  ENT.close();
}
