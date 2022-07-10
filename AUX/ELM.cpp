#include "ELM.h"

CEl::CEl(const int &_NEl, const int &_DegreeApproxPol, const int &_NMa, tvMa &vMa, const int &_NNo, const V1I &vNNo, tvNo &vNo, const int &_NQuad, tvQuadratura &_Quad, CCo &_Co)
{
  NEl = _NEl;
  dimension = 2;
  DegreeApproxPol = _DegreeApproxPol;
  NNo = _NNo;
  dDegreeApproxPol = DegreeApproxPol - 1;
  NNEl1m = c_NNEl1m();

  for (int a = 0; a < this->NNo; a++)
  {
    No.push_back(&vNo[vNNo[a]]);
    no.push_back(vNo[vNNo[a]].r_NumberNo());
    No[a]->w_FP(1);
    No[a]->I_Elm(NEl, a);
  }

  NMa = _NMa;
  Ma = &vMa[NMa]; // talvez troca por _NMa ..
  NQuad = _NQuad;
  Quad = &_Quad[NQuad];
  Config = &_Co;
  Xsi.resize(NNo);
  no_E.resize(NNo);
  no_S.resize(NNo);
  no_Sig.resize(NNo);
  no_ue.resize(NNo, 0.0);

  for (int i = 0; i < NNo; ++i)
  {
    Xsi[i].resize(dimension, 0.0);
    no_E[i].resize(dimension * dimension, 0.0);
    no_S[i].resize(dimension * dimension, 0.0);
    no_Sig[i].resize(dimension * dimension, 0.0);
  }

  Fint.resize(dimension * NNo, 0.0);
  Hessiana.resize(dimension * NNo, V1D(dimension * NNo, 0.0));
  Mass.resize(dimension * NNo, V1D(dimension * NNo, 0.0));
  MRP.resize(NNo, V1D(NNo, 0.0));
  PHI.resize(r_Quad()->r_NumberQuad());
  DPHI.resize(r_Quad()->r_NumberQuad());
  X.resize(r_Quad()->r_NumberQuad());
  A0.resize(r_Quad()->r_NumberQuad());
  J0.resize(r_Quad()->r_NumberQuad(), 0.0);
  invA0.resize(r_Quad()->r_NumberQuad());
  Y.resize(r_Quad()->r_NumberQuad());
  A1.resize(r_Quad()->r_NumberQuad());
  J1.resize(r_Quad()->r_NumberQuad(), 0.0);
  A.resize(r_Quad()->r_NumberQuad());
  J.resize(r_Quad()->r_NumberQuad(), 0.0);
  Sig.resize(r_Quad()->r_NumberQuad());
  C.resize(r_Quad()->r_NumberQuad());
  E.resize(r_Quad()->r_NumberQuad());
  S.resize(r_Quad()->r_NumberQuad());
  CC.resize(r_Quad()->r_NumberQuad());
  ue.resize(r_Quad()->r_NumberQuad(), 0.0);
  fint.resize(r_Quad()->r_NumberQuad());
  h.resize(r_Quad()->r_NumberQuad());
  m.resize(r_Quad()->r_NumberQuad());
  mrp.resize(r_Quad()->r_NumberQuad());
  vrp.resize(r_Quad()->r_NumberQuad());

  for (int i = 0; i < r_Quad()->r_NumberQuad(); ++i)
  {
    PHI[i].resize(NNo, 0.0);
    DPHI[i].resize(NNo, V1D(dimension, 0.0));
    X[i].resize(dimension, 0.0);
    A0[i].resize(dimension, V1D(dimension, 0.0));
    invA0[i].resize(dimension, V1D(dimension, 0.0));
    Y[i].resize(dimension, 0.0);
    A1[i].resize(dimension, V1D(dimension, 0.0));
    A[i].resize(dimension, V1D(dimension, 0.0));
    Sig[i].resize(dimension * dimension, 0.0);
    C[i].resize(dimension, V1D(dimension, 0.0));
    E[i].resize(dimension * dimension, 0.0);
    S[i].resize(dimension * dimension, 0.0);
    CC[i].resize(dimension * dimension, V1D(dimension * dimension, 0.0));
    fint[i].resize(dimension * NNo, 0.0);
    h[i].resize(dimension * NNo, V1D(dimension * NNo, 0.0));
    m[i].resize(dimension * NNo, V1D(dimension * NNo, 0.0));
    mrp[i].resize(NNo, V1D(NNo, 0.0));
    vrp[i].resize(NNo, 0.0);
  }

  Xsi = c_Xsi(NNo, DegreeApproxPol, dimension);

  //================== VERIFICAÇÃO =======================
  /*
  std::cout << "Número de pontos - Coordenadas homogêneas -  ξ1  ξ2" << "\n";
  for (int i = 0; i < NNo; i++)
  {
  std::cout << "i=" << i << " Xsi0=" << Xsi[i][0] << " Xsi1=" << Xsi[i][1] << std::endl;
  }
  */
  // =====================================================

  Ue = 0.0;
  for (int i = 0; i < r_Quad()->r_NumberQuad(); i++)
  {
    //std::cout << " Nó " << i << " da Quadratura" << "\n";
    PHI[i] = calc_PHI(DegreeApproxPol, NNo, dimension, r_Quad()->r_Xsi(i));
    DPHI[i] = calc_DPHI(DegreeApproxPol, NNo, dimension, r_Quad()->r_Xsi(i));
    A0[i] = calc_An(DPHI[i], No, 0); //std::cout << "A0: [" << A0[i][0][0]  << ", " << A0[i][0][1]  <<", " << A0[i][1][0] << ", "<<  A0[i][1][1]  <<"\n";
    J0[i] = Det(A0[i]);  //std:: cout << "J0: " << J0[i] << "\n";
    X[i] = calc_Xn(dimension, PHI[i], No, 0);
    invA0[i] = invEIGEN(A0[i]); //std::cout << "invA0: [" << invA0[i][0][0]  << ", " << invA0[i][0][1]  <<", " << invA0[i][1][0] << ", "<<  invA0[i][1][1]  <<"\n";
    //std::cout << "=============================" << "\n";
  }
}

void CEl::ite()
{
  // DINAMICO -> ACELERAÇÃO MEDIA CTE (Beta=0.25 e Gamma =0.5)
  // DINAMICO -> ACELERAÇÃO LINEAR    (Beta=1/6 e Gamma =0.5)
  double beta_ = 1.0;
  double deltat_ = 1.0;
  for (int i = 0; i < dimension * NNo; i++)
  {
    Fint[i] = 0.0;
    for (int j = 0; j < dimension * NNo; j++)
    {
      Hessiana[i][j] = 0.0;
      Mass[i][j]=0.0; // ====== DINAMICO ======
    }
  }
  Ue = 0.0;
  for (int i = 0; i < r_Quad()->r_NumberQuad(); i++)
  {
    Y[i] = calc_Xn(dimension, PHI[i], No, 1);
    A1[i] = calc_An(DPHI[i], No, 1); // Gradiente de tranformação final
    J1[i] = Det(A1[i]);
    A[i] = MM(A1[i], invA0[i]);
    J[i]= Det(A[i]);                          //std::cout << "J:" << J[i] << "\n"; 
    C[i] = MtM(A[i], A[i]);
    E[i] = calc_Ei(C[i]);                                                                             // executa a cada iteracao para cada ponto de gauss e armazena Si=due_dE, e CCi=d2ue_dEdE
    CC[i] = calc_CC(dimension, Ma);                                                                   // Calculo do tensor contitutivo de Piola-Kirchhoff de segunda espécie;
    S[i] = calc_S(dimension, E[i], CC[i]);
    Sig[i]= calc_Sig(dimension, A[i], J[i], S[i]);
    calc_Fint_H_M(PHI[i], DPHI[i], invA0[i], A[i], S[i], CC[i], fint[i], h[i], m[i], mrp[i], vrp[i], A1[i]); // Calculo do vetor de forças internas e da matriz Hesiana
    ue[i] = c_ue(E[i], Ma);
    Ue = Ue + ue[i] * J0[i] * r_Quad()->r_W(i);
    for (int a = 0; a < NNo; a++)
    {
      for (int z = 0; z < dimension; z++)
      {
        Fint[dimension * a + z] += fint[i][dimension * a + z] * J0[i] * r_Quad()->r_W(i);
        for (int l = 0; l < NNo; l++)
        {
          for (int k = 0; k < dimension; k++)
          {
            Hessiana[dimension * a + z][dimension * l + k] += h[i][dimension * a + z][dimension * l + k] * J0[i] * r_Quad()->r_W(i);
            // for(int linha=0; linha <NNo; linha++)
            // {
            //   for (int coluna = 0; coluna< NNo; coluna++)
            //   {
            //   std::cout << Hessiana[linha][coluna];
            //   }
            // }
            if (z == k)
            Mass[dimension * a + z][dimension * l + k] += (m[i][dimension * a + z][dimension * l + k] / (beta_ * deltat_ * deltat_)) * J0[i] * r_Quad()->r_W(i); /**/ // mass contribution
          }
        }
      }
    }
    for (int a = 0; a < NNo; a++)
    {
      // VRP[a]+=VRP[i][a]*J0[i]*r_Quad()->r_W(i);
      for (int l = 0; l < NNo; l++)
      {
        MRP[a][l] += (mrp[i][a][l]) * J0[i] * r_Quad()->r_W(i);
      }
    }
    if (J[i]<0)
    {
      std::cout << "Jacobiano invertido: " << J[i] <<"\n";
    }
  }



  
  // metodo dos residuos ponderados (Rodolfo)

  // metodo minimos quadrados (Jeferson)

  // metdo minimosquadrados + 3 pontos (Coda)

  // metodo nodal WWW (NO A NO + MEDIA NO CONTORNO)

  for (int l = 0; l < NNo; l++)
  {
    cpto(DegreeApproxPol, NNo, dimension, Xsi[l], No, Ma, no_E[l], no_S[l], no_Sig[l], no_ue[l]);
  }

  // == ALTERAR PARA UMA VARIAVEL DE LEITURA E/OU PASSANDO COMO ARGUMENTO NA FUNÇÃO
    int impRel;
    impRel = 1;
    if (impRel == 0)
    {
      Imp_Rel_Geral();
    }
}

void cpto(const int &DegreeApproxPol, const int &NNo, const int &dimension, const V1D XsiNo, std::vector<CNo *> No, CMa *Ma, V1D &E, V1D &S, V1D &Sig, double &ue)
{
  V1D PHI;   ///< Valor das funções de forma nos pontos de integração
  V2D DPHI;  ///< Valor da primeira derivada das funcoes de forma nos pontos de integração
  V1D X;     ///< Posição do ponto de integração na configuração inicial
  V2D A0;    ///< Gradiente de tranformação inicial
  double J0; ///< Determinante do Gradiente de tranformação inicial
  V2D invA0; ///< Matriz inversa da Gradiente de tranformação inicial
  V1D Y;     ///< Posição do ponto de integração na configuração final
  V2D A1;    ///< Gradiente de tranformação final
  double J1; ///< Determinante do Gradiente de tranformação final
  double J;  ///< Determinante do Gradiente de tranformação
  V2D A;     ///< Gradiente de transformação total
  V2D C;     ///< Tensor de alongamento ou estiramento a direita de Cauchy-Green
  V2D CC;    ///< Tensor constitutivo de Piola Kirchhoff de segunda espécie
  V1D fint;  ///< Vetor de forças específica* internas
  V2D h;     ///< Matrix Hesiana específica*
  V2D m;     ///< Matrix de massa específica*
  V2D mrp;   ///< MRP <Matriz aproximadora por residuos ponderados local>
  V1D vrp;   ///< VRP <Matriz aproximadora por residuos ponderados local>

  PHI.resize(NNo, 0.0);
  DPHI.resize(NNo, V1D(dimension, 0.0));
  X.resize(dimension, 0.0);
  A0.resize(dimension, V1D(dimension, 0.0));
  J0 = 0.0;
  invA0.resize(dimension, V1D(dimension, 0.0));
  Y.resize(dimension, 0.0);
  A1.resize(dimension, V1D(dimension, 0.0));
  J1 = 0.0;
  J  = 0.0;
  A.resize(dimension, V1D(dimension, 0.0));
  C.resize(dimension, V1D(dimension, 0.0));
  CC.resize(dimension * dimension, V1D(dimension * dimension, 0.0));
  fint.resize(dimension * NNo, 0.0);
  h.resize(dimension * NNo, V1D(dimension * NNo, 0.0));
  m.resize(dimension * NNo, V1D(dimension * NNo, 0.0));
  mrp.resize(NNo, V1D(NNo, 0.0));
  vrp.resize(NNo, 0.0);
  // std::cout << " " << m[0][0] << std::endl;
  //   std::cout << "Xsi " << XsiNo[0] << " " << XsiNo[1] << " " << XsiNo[2] << std::endl;
  //   std::cout << "FF"<< std::endl; for (int jj=0; jj<NN; jj++) {   for (int ii=0; ii<NN; ii++) { std::cout << FF[jj][ii] << " "; } std::cout<<std::endl; }
  PHI = calc_PHI(DegreeApproxPol, NNo, dimension, XsiNo);   // std::cout << "PHI"<< std::endl; for (int ii=0; ii<PHI.size(); ii++) { std::cout << PHI[0] << " "; } std::cout<<std::endl;
  DPHI = calc_DPHI(DegreeApproxPol, NNo, dimension, XsiNo); // std::cout << "DPHI"<< std::endl; for (int jj=0; jj<d; jj++) {   for (int ii=0; ii<DPHI.size(); ii++) { std::cout << DPHI[ii][jj] << " "; } std::cout<<std::endl; }
  A0 = calc_An(DPHI, No, 0);
  J0 = Det(A0);
  X = calc_Xn(dimension, PHI, No, 0); //  std::cout << "X " << X[0] << " " << X[1] << " " << X[2] << std::endl;
  invA0 = invEIGEN(A0);
  Y = calc_Xn(dimension, PHI, No, 1); //  std::cout << "Y " << Y[0] << " " << Y[1] << " " << Y[2] << std::endl;
  A1 = calc_An(DPHI, No, 1);          // Gradiente de tranformação final
  J1 = Det(A1);
  A = MM(A1, invA0);                                               // std::cout << "A " << A[0][0] << " " << A[0][1] << " " << A[0][2] << " " << A[1][0] << " " << A[1][1] << " " << A[1][2] << " " << A[2][0] << " " << A[2][1] << " " << A[2][2] << std::endl;
  J = Det(A);                          //std::cout << "J: "<< J << "\n"; 
  C = MtM(A, A);                                                   // std::cout << "C " << C[0][0] << " " << C[0][1] << " " << C[0][2] << " " << C[1][0] << " " << C[1][1] << " " << C[1][2] << " " << C[2][0] << " " << C[2][1] << " " << C[2][2] << std::endl;
  E = calc_Ei(C);                                                  // std::cout << "E " << E[0][0] << " " << E[0][1] << " " << E[0][2] << " " << E[1][0] << " " << E[1][1] << " " << E[1][2] << " " << E[2][0] << " " << E[2][1] << " " << E[2][2] << std::endl;
  CC = calc_CC(dimension, Ma);                                     // Calculo do tensor contitutivo de Piola-Kirchhoff de segunda espécie;
  S = calc_S(dimension, E, CC);                                    // Calculo da Tensões de Green;
  Sig = calc_Sig(dimension, A, J, S);  
  calc_Fint_H_M(PHI, DPHI, invA0, A, S, CC, fint, h, m, mrp, vrp, A1); // Calculo do vetor de forças internas e da matriz Hesiana
  ue = c_ue(E, Ma);
}

CEl::~CEl()
{
  No.clear();
  no.clear();
  Xsi.clear();
  Fint.clear();
  Hessiana.clear();
  PHI.clear();
  DPHI.clear();
  X.clear();
  A0.clear();
  J0.clear();
  invA0.clear();
  Y.clear();
  A1.clear();
  J1.clear();
  A.clear();
  C.clear();
  E.clear();
  CC.clear();
  S.clear();
  ue.clear();
  fint.clear();
  h.clear();
  m.clear();
}

int CEl::c_NNEl1m()
{
  int nn1m;
  nn1m = (dDegreeApproxPol + 1) * ((dDegreeApproxPol + 1) + 1) / 2; // numero de nos do elemento triangular com grau do polinomio uma vez menor
  return nn1m;
}

V2D c_Xsi(const int &NNo, const int &DegreeApproxPol, const int &dimension)
{
  V2D Xsi;
  Xsi.resize(NNo, V1D(dimension, 0.0));
  // NN que é dado pode também ser calculado por: NN=(Pe+1)*(2+Pe)/2;
  // Calcula as Coordenadas Adimensionais dos Nós
  Xsi[0][0] = 0.0;
  Xsi[0][1] = 0.0;
  Xsi[1][0] = 1.0;
  Xsi[1][1] = 0.0;
  Xsi[2][0] = 0.0;
  Xsi[2][1] = 1.0;
  Xsi[3][0] = 1.0 / 3.0;
  Xsi[3][1] = 0.0;
  Xsi[4][0] = 2.0 / 3.0;
  Xsi[4][1] = 0.0;
  Xsi[5][0] = 2.0 / 3.0;
  Xsi[5][1] = 1.0 / 3.0;
  Xsi[6][0] = 1.0 / 3.0;
  Xsi[6][1] = 2.0 / 3.0;
  Xsi[7][0] = 0.0;
  Xsi[7][1] = 2.0 / 3.0;
  Xsi[8][0] = 0.0;
  Xsi[8][1] = 1.0 / 3.0;
  Xsi[9][0] = 1.0 / 3.0;
  Xsi[9][1] = 1.0 / 3.0;

  return Xsi;
}

V1D calc_PHI(const int &DegreeApproxPol, const int &NNo, const int dimension, const V1D &Xsi)
{
  V1D PHIi;
  PHIi.resize(NNo, 0.0);

  switch (DegreeApproxPol)
  {
  case 3:
  {
    PHIi[0] = 1 - (11 * Xsi[1]) / 2. + 9 * pow(Xsi[1], 2) - (9 * pow(Xsi[1], 3)) / 2. - (11 * Xsi[0]) / 2. + 18 * Xsi[1] * Xsi[0] - (27 * pow(Xsi[1], 2) * Xsi[0]) / 2. + 9 * pow(Xsi[0], 2) - (27 * Xsi[1] * pow(Xsi[0], 2)) / 2. - (9 * pow(Xsi[0], 3)) / 2.;
    PHIi[1] = Xsi[0] - (9 * pow(Xsi[0], 2)) / 2. + (9 * pow(Xsi[0], 3)) / 2.;
    PHIi[2] = Xsi[1] - (9 * pow(Xsi[1], 2)) / 2. + (9 * pow(Xsi[1], 3)) / 2.;
    PHIi[3] = 9 * Xsi[0] - (45 * Xsi[1] * Xsi[0]) / 2. + (27 * pow(Xsi[1], 2) * Xsi[0]) / 2. - (45 * pow(Xsi[0], 2)) / 2. + 27 * Xsi[1] * pow(Xsi[0], 2) + (27 * pow(Xsi[0], 3)) / 2.;
    PHIi[4] = (-9 * Xsi[0]) / 2. + (9 * Xsi[1] * Xsi[0]) / 2. + 18 * pow(Xsi[0], 2) - (27 * Xsi[1] * pow(Xsi[0], 2)) / 2. - (27 * pow(Xsi[0], 3)) / 2.;
    PHIi[5] = (-9 * Xsi[1] * Xsi[0]) / 2. + (27 * Xsi[1] * pow(Xsi[0], 2)) / 2.;    
    PHIi[6] = (-9 * Xsi[1] * Xsi[0]) / 2. + (27 * pow(Xsi[1], 2) * Xsi[0]) / 2.;
    PHIi[7] = (-9 * Xsi[1]) / 2. + 18 * pow(Xsi[1], 2) - (27 * pow(Xsi[1], 3)) / 2. + (9 * Xsi[1] * Xsi[0]) / 2. - (27 * pow(Xsi[1], 2) * Xsi[0]) / 2.;
    PHIi[8] = 9 * Xsi[1] - (45 * pow(Xsi[1], 2)) / 2. + (27 * pow(Xsi[1], 3)) / 2. - (45 * Xsi[1] * Xsi[0]) / 2. + 27 * pow(Xsi[1], 2) * Xsi[0] + (27 * Xsi[1] * pow(Xsi[0], 2)) / 2.;
    PHIi[9] = 27 * Xsi[1] * Xsi[0] - 27 * pow(Xsi[1], 2) * Xsi[0] - 27 * Xsi[1] * pow(Xsi[0], 2);
    break;
  }
  default:
  {
    std::cout << "Aproximação diferente, corrigir para aproximação cúbica (3)!."
              << "\n";
    break;
  }
  break;
  }
  return PHIi;
}

V2D calc_DPHI(const int &DegreeApproxPol, const int &NNo, const int dimension, const V1D &Xsi)
{
  V2D DPHIi;
  DPHIi.resize(NNo, V1D(dimension, 0.0));

  switch (DegreeApproxPol)
  {
  case 3:
  {
    DPHIi[0][0] = -5.5 + 18 * Xsi[1] - (27 * pow(Xsi[1], 2)) / 2. + 18 * Xsi[0] - 27 * Xsi[1] * Xsi[0] - (27 * pow(Xsi[0], 2)) / 2.;
    DPHIi[1][0] = 1 - 9 * Xsi[0] + (27 * pow(Xsi[0], 2)) / 2.;
    DPHIi[2][0] = 0;
    DPHIi[3][0] = 9 - (45 * Xsi[1]) / 2. + (27 * pow(Xsi[1], 2)) / 2. - 45 * Xsi[0] + 54 * Xsi[1] * Xsi[0] + (81 * pow(Xsi[0], 2)) / 2.;
    DPHIi[4][0] = -4.5 + (9 * Xsi[1]) / 2. + 36 * Xsi[0] - 27 * Xsi[1] * Xsi[0] - (81 * pow(Xsi[0], 2)) / 2.;
    DPHIi[5][0] = (-9 * Xsi[1]) / 2. + 27 * Xsi[1] * Xsi[0];
    DPHIi[6][0] = (-9 * Xsi[1]) / 2. + (27 * pow(Xsi[1], 2)) / 2.;
    DPHIi[7][0] = (9 * Xsi[1]) / 2. - (27 * pow(Xsi[1], 2)) / 2.;
    DPHIi[8][0] = (-45 * Xsi[1]) / 2. + 27 * pow(Xsi[1], 2) + 27 * Xsi[1] * Xsi[0];
    DPHIi[9][0] = 27 * Xsi[1] - 27 * pow(Xsi[1], 2) - 54 * Xsi[1] * Xsi[0];

    DPHIi[0][1] = -5.5 + 18 * Xsi[1] - (27 * pow(Xsi[1], 2)) / 2. + 18 * Xsi[0] - 27 * Xsi[1] * Xsi[0] - (27 * pow(Xsi[0], 2)) / 2.;
    DPHIi[1][1] = 0;
    DPHIi[2][1] = 1 - 9 * Xsi[1] + (27 * pow(Xsi[1], 2)) / 2.;
    DPHIi[3][1] = (-45 * Xsi[0]) / 2. + 27 * Xsi[1] * Xsi[0] + 27 * pow(Xsi[0], 2);
    DPHIi[4][1] = (9 * Xsi[0]) / 2. - (27 * pow(Xsi[0], 2)) / 2.;
    DPHIi[5][1] = (-9 * Xsi[0]) / 2. + (27 * pow(Xsi[0], 2)) / 2.;
    DPHIi[6][1] = (-9 * Xsi[0]) / 2. + 27 * Xsi[1] * Xsi[0]; 
    DPHIi[7][1] = -4.5 + 36 * Xsi[1] - (81 * pow(Xsi[1], 2)) / 2. + (9 * Xsi[0]) / 2. - 27 * Xsi[1] * Xsi[0];  
    DPHIi[8][1] = 9 - 45 * Xsi[1] + (81 * pow(Xsi[1], 2)) / 2. - (45 * Xsi[0]) / 2. + 54 * Xsi[1] * Xsi[0] + (27 * pow(Xsi[0], 2)) / 2.;
    DPHIi[9][1] = 27 * Xsi[0] - 54 * Xsi[1] * Xsi[0] - 27 * pow(Xsi[0], 2);
  
    break;
  }
  default:
  {
    std::cout << "Aproximação diferente, corrigir para aproximação cúbica (3)!."
              << "\n";
    break;
  }
  break;
  }
  return DPHIi;
}

V1D calc_Xn(const int &dimension, const V1D &PHI, std::vector<CNo *> No, const int estado)
{
  V1D xi;
  xi.resize(dimension, 0.0);
  int NN = PHI.size();
  for (int m = 0; m < dimension; m++)
  {
    xi[m] = 0.0;
    for (int L = 0; L < NN; L++)
    {
      xi[m] += PHI[L] * No[L]->r_X(estado, m);
    }
  }
  return xi;
}

V2D calc_An(const V2D &DPHI, std::vector<CNo *> No, const int estado)
{
  int NN = DPHI.size();
  int d = DPHI[0].size();
  V2D Ani;
  Ani.resize(d, V1D(d, 0.0));
  for (int m = 0; m < d; m++)
  {
    for (int n = 0; n < d; n++)
    {
      Ani[m][n] = 0.0;
      for (int L = 0; L < DPHI.size(); L++)
      {
        Ani[m][n] += DPHI[L][n] * No[L]->r_X(estado, m);
      }
    }
  }
  return Ani;
}

V1D calc_Ei(const V2D &Ci)
{
  V1D Ei;
  Ei.resize(Ci.size() * Ci[0].size(), 0.0);
  Ei[0] = (Ci[0][0] - 1.0) / 2.0;
  Ei[1] = (Ci[1][1] - 1.0) / 2.0;
  //Ei[2] = (Ci[0][1]);
  //Ei[3] = (Ci[1][0]);
  Ei[2] = (Ci[0][1]) / 2.0;
  Ei[3] = (Ci[1][0]) / 2.0;
  return Ei;
}

V2D calc_CC(const int &dimension, CMa *Ma)
{
  V2D CCi;
  CCi.resize(dimension * dimension, V1D(dimension * dimension, 0.0));
  //EPD
  CCi[0][0] = (2 * Ma->r_LameM() + Ma->r_LameL());
  CCi[0][1] = Ma->r_LameL();
  CCi[1][0] = Ma->r_LameL();
  CCi[1][1] = (2 * Ma->r_LameM() + Ma->r_LameL());
  /// verificar direito estes termos:
  //CCi[2][2] = (Ma->r_LameM())*2;
  //CCi[3][3] = (Ma->r_LameM())*2;
  CCi[2][2] = Ma->r_LameM();
  CCi[2][3] = Ma->r_LameM();   
  CCi[3][2] = Ma->r_LameM();   
  CCi[3][3] = Ma->r_LameM();

  //EPT
  //CCi[0][0] = (Ma->r_E0()/(1-pow(Ma->r_Nu0(),2)));
  //CCi[0][1] = (Ma->r_E0()/(1-pow(Ma->r_Nu0(),2)))*Ma->r_Nu0();
  //CCi[1][0] = (Ma->r_E0()/(1-pow(Ma->r_Nu0(),2)))*Ma->r_Nu0();
  //CCi[1][1] = (Ma->r_E0()/(1-pow(Ma->r_Nu0(),2)));
  // verificar direito estes termos:
  //CCi[2][2] = Ma->r_Et();
  //CCi[2][3] = Ma->r_Et();
  //CCi[3][2] = Ma->r_Et();
  //CCi[3][3] = Ma->r_Et();

  return CCi;
}
  
V1D calc_S(const int &dimension, const V1D &Ei, V2D &CCi)
{
  V1D Si;
  Si.resize(dimension * dimension, 0.0);
  Si[0] = CCi[0][0] * Ei[0] + CCi[0][1] * Ei[1];
  Si[1] = CCi[1][0] * Ei[0] + CCi[1][1] * Ei[1];
  Si[2] = CCi[2][2] * Ei[2] + CCi[2][3] * Ei[3]; 
  Si[3] = CCi[3][2] * Ei[2] + CCi[3][3] * Ei[3]; 
  return Si;
}

double c_ue(const V1D &E, CMa *Ma)
{
  double ue;
  ue = (1.0 / 2.0) *
       ((2.0 * Ma->r_LameM() + Ma->r_LameL()) * (E[0] * E[0] + E[1] * E[1]) +
        (Ma->r_LameL()) * (E[0] * E[1]) +
        (Ma->r_LameM()) * ((E[2] + E[3]) * (E[2] + E[3]))); //EPD
  return ue;
}

void calc_Fint_H_M(const V1D &PHI, const V2D &DPHI, const V2D &Di, const V2D &Ai, const V1D &Si, const V2D &CCi, V1D &finti, V2D &hi, V2D &mi, V2D &mrpi, V1D &vrpi, V2D const A1)
{
  int NN = DPHI.size();
  int d = DPHI[0].size();
  V4D dA1_dY;
  dA1_dY.resize(NN, V3D(d, V2D(d, V1D(d, 0.0))));
  V4D dA_dY;
  dA_dY.resize(NN, V3D(d, V2D(d, V1D(d, 0.0))));
  V1D dE_dYaz;
  dE_dYaz.resize(d * d, 0.0);
  V1D dE_dYlk;
  dE_dYlk.resize(d * d, 0.0);
  V1D d2E_dYdY;
  d2E_dYdY.resize(d * d, 0.0);
  V1D CCi_dE_dYlk;
  CCi_dE_dYlk.resize(d, 0.0);
  V2D AdA_dY;
  AdA_dY.resize(d, V1D(d, 0.0));
  V2D dA_dY_dA_dY_1;
  dA_dY_dA_dY_1.resize(d, V1D(d, 0.0));
  V2D dA_dY_dA_dY_2;
  dA_dY_dA_dY_2.resize(d, V1D(d, 0.0));

  //> DINAMICO A SER IMPLEMENTADO
  double density = 1.0;
  //=============================

  for (int L = 0; L < NN; L++)
  {
    for (int m = 0; m < d; m++)
    {
      for (int n = 0; n < d; n++)
      {
        dA1_dY[L][m][m][n] = DPHI[L][n];
      }
    }
  }
  for (int a = 0; a < NN; a++)
  {
    for (int z = 0; z < d; z++)
    {
      dA_dY[a][z] = MM(dA1_dY[a][z], Di);
    }
  }

  // Calcula Matriz hessiana e fint no ponto de integração
  for (int a = 0; a < NN; a++)
  {
    for (int z = 0; z < d; z++)
    {
      //AdA_dY = MtM(Ai, dA_dY[a][z]);
      AdA_dY = MtM(dA_dY[a][z],Ai);
      dE_dYaz[0] = (1.0 / 2.0) * (AdA_dY[0][0] + AdA_dY[0][0]);
      dE_dYaz[1] = (1.0 / 2.0) * (AdA_dY[1][1] + AdA_dY[1][1]);
      dE_dYaz[2] = (1.0 / 2.0) * (AdA_dY[1][0] + AdA_dY[0][1]);
      dE_dYaz[3] = (1.0 / 2.0) * (AdA_dY[0][1] + AdA_dY[1][0]);
      finti[2 * a + z] = Si[0] * dE_dYaz[0] + Si[1] * dE_dYaz[1] + Si[2] * dE_dYaz[2] + Si[3] * dE_dYaz[3];
      for (int l = 0; l < NN; l++)
      {
        for (int k = 0; k < d; k++)
        {
          //AdA_dY = MtM(Ai, dA_dY[l][k]);
          AdA_dY = MtM(dA_dY[l][k],Ai);
          dE_dYlk[0] = (1.0 / 2.0) * (AdA_dY[0][0] + AdA_dY[0][0]);
          dE_dYlk[1] = (1.0 / 2.0) * (AdA_dY[1][1] + AdA_dY[1][1]);
          dE_dYlk[2] = (1.0 / 2.0) * (AdA_dY[1][0] + AdA_dY[0][1]);
          dE_dYlk[3] = (1.0 / 2.0) * (AdA_dY[0][1] + AdA_dY[1][0]);
          dA_dY_dA_dY_1 = MtM(dA_dY[a][z], dA_dY[l][k]);
          dA_dY_dA_dY_2 = MtM(dA_dY[l][k], dA_dY[a][z]);
          d2E_dYdY[0] = (1.0 / 2.0) * (dA_dY_dA_dY_1[0][0] + dA_dY_dA_dY_2[0][0]);
          d2E_dYdY[1] = (1.0 / 2.0) * (dA_dY_dA_dY_1[1][1] + dA_dY_dA_dY_2[1][1]);
          d2E_dYdY[2] = (1.0 / 2.0) * (dA_dY_dA_dY_1[0][1] + dA_dY_dA_dY_2[0][1]);
          d2E_dYdY[3] = (1.0 / 2.0) * (dA_dY_dA_dY_1[1][0] + dA_dY_dA_dY_2[1][0]);
          CCi_dE_dYlk = MV(CCi, dE_dYlk);
          hi[2 * a + z][2 * l + k] = dE_dYaz[0] * CCi_dE_dYlk[0] +
                                     dE_dYaz[1] * CCi_dE_dYlk[1] +
                                     dE_dYaz[2] * CCi_dE_dYlk[2] +
                                     dE_dYaz[3] * CCi_dE_dYlk[3] +
                                     (Si[0] * d2E_dYdY[0] +
                                      Si[1] * d2E_dYdY[1] +
                                      Si[2] * d2E_dYdY[2] +
                                      Si[3] * d2E_dYdY[3]);
          if (z == k)
            mi[2 * a + z][2 * l + k] = density * PHI[a] * PHI[l]; // mass contribution
        }
      }
    }
  }
  for (int a = 0; a < NN; a++)
  {
    for (int l = 0; l < NN; l++)
    {
      mrpi[a][l] = PHI[a] * PHI[l]; // matriz de interpolacao via residuos ponderados
    }
    vrpi[a] = PHI[a]; // vetor de interpolacao via residuos ponderados
  }
  //std::cout << "A :[" << Ai[0][0] << ", "<< Ai[1][1]  <<", "<< Ai[0][1] <<", "<< Ai[1][0] << "]\n";
}

void Le_Elementos(tvEl &El, tvNo &No, tvMa &Ma, tvQuadratura &Qua, CCo &Co)
{
  int i, nEl;
  char s[1000];
  std::ifstream Ent; // para leitura
  std::string NAr = Co.r_ArqELM();
  Ent.open(NAr.c_str());
  Ent.getline(s, 1000);
  Ent >> nEl;
  Ent.getline(s, 1000);
  Ent.getline(s, 1000);
  for (i = 0; i < nEl; i++)
  {
    V1I vNNo;
    int nno = 10; ///< Numero de nos do elemento;
    int D = 2;    ///< Dimensão do elemento: 1D, 2D ou 3D;
    int Pe = 3;   ///< Aproximacao fisica do elemento: 1=linear 2=quadratico ...;
    int ma = 0;   ///< Indice do material que compoe o elemento;
    int qua = 0;  ///< Indice da quadratura a ser usada para integrar o elemento;
    int I;        ///< Indice do elemento;
    int Type_Element;
    Ent >> Type_Element;
    vNNo.resize(nno, 0);
    for (int j = 0; j < nno; j++)
    {
      Ent >> vNNo[j];
    }
    Ent >> I;
    Ent.getline(s, 1000);
    CEl LeEl(i, Pe, ma, Ma, nno, vNNo, No, qua, Qua, Co);
    El.push_back(LeEl);
  }
  Ent.close();
}

V1D calc_Sig (const int& dimension, const V2D& A,const int& J,const V1D& S)
{
  V2D mSig;
  V2D mS;
  V1D Sig;
  V2D At;
  V2D mSxAt;
  double Jacc;
  mSig.resize(dimension, V1D(dimension,0.0));
  mSxAt.resize(dimension, V1D(dimension,0.0));
  At.resize(dimension, V1D(dimension,0.0));
  mS.resize(dimension, V1D(dimension,0.0));
  Sig.resize(dimension*dimension, 0.0);
  mS[0][0]=S[0];
  mS[1][1]=S[1];
  mS[0][1]=S[2];
  mS[1][0]=S[3];
  At=Mt(A);         //std::cout << " A: " << "[" << A[0][0] << ", " << A[0][1] << "][" << A[1][0] <<", " << A[1][1] << "]" <<"\n";
  mSxAt=MM(mS, Mt(A));
  mSig=MM(A,mSxAt);
  Jacc=Det(A);
  Sig[0]=mSig[0][0]/Jacc;
  Sig[1]=mSig[1][1]/Jacc;
  Sig[2]=mSig[0][1]/Jacc;
  Sig[3]=mSig[1][0]/Jacc; //std::cout << "Jacc: " << Jacc << "\t [" << Sig[0] << ", " << Sig[1] << ", " << Sig[2] <<", " << Sig[3] << "]" <<"\n"; std::cout << "-----------------------------------------------------" << "\n";
  return Sig;
}

void CEl::Imp_Rel_Geral()
{
  std::string NArq1;
  NArq1.append(Config->r_ArqREL());
  std::fstream Sai1;
  Sai1.open(NArq1.c_str(), std::fstream::in | std::fstream::out | std::fstream::app); ///< Dados Calculados
  Sai1 << "Iteração ="
       << "X" << std::endl;
  Sai1 << "Elemento =" << NEl << std::endl;
  Sai1 << "Grau do Polinomio aproximador =" << DegreeApproxPol << std::endl;
  Sai1 << "Numero de nós =" << NNo << std::endl;
  Sai1 << "Coordenadas Adimensionais dos Nós" << std::endl;
  for (int I = 0; I < Xsi.size(); I++)
  {
    Sai1 << I << " ";
    for (int J = 0; J < Xsi[0].size(); J++)
    {
      Sai1 << "Xsi[" << I << "][" << J << "]=" << Xsi[I][J] << "   ";
    }
    Sai1 << std::endl;
  }
  for (int i = 0; i < r_Quad()->r_NumberQuad(); i++)
  {
    Sai1 << "Ponto de integracao No. =" << i << std::endl;
    Sai1 << "Valores das funcoes de forma no ponto de Integracao atual(PHI) :" << std::endl;
    for (int L = 0; L < PHI[0].size(); L++)
    {
      Sai1 << PHI[i][L] << " ";
    }
    Sai1 << std::endl;
    Sai1 << "Valores das derivadas direcionais das funcoes de forma no ponto de Integracao atual(DPHI):" << std::endl;
    for (int M = 0; M < dimension; M++)
    {
      for (int L = 0; L < DPHI[0].size(); L++)
      {
        Sai1 << DPHI[i][L][M] << " ";
      }
      Sai1 << std::endl;
    }
    Sai1 << " X= " << std::endl;
    for (int m = 0; m < dimension; m++)
    {
      Sai1 << X[i][m] << " ";
    }
    Sai1 << std::endl;
    Sai1 << " J0= " << J0[i] << std::endl;
    Sai1 << " Y= " << std::endl;
    for (int m = 0; m < dimension; m++)
    {
      Sai1 << Y[i][m] << " ";
    }
    Sai1 << std::endl;
    Sai1 << "Matriz A1= " << std::endl;
    for (int m = 0; m < dimension; m++)
    {
      for (int n = 0; n < dimension; n++)
      {
        Sai1 << A1[i][m][n] << " ";
      }
      Sai1 << std::endl;
    }
    Sai1 << " J1= " << J1[i] << std::endl;
    Sai1 << " A= " << std::endl;
    for (int m = 0; m < dimension; m++)
    {
      for (int n = 0; n < dimension; n++)
      {
        Sai1 << A[i][m][n] << " ";
      }
      Sai1 << std::endl;
    }
    Sai1 << " C= " << std::endl;
    for (int m = 0; m < dimension; m++)
    {
      for (int n = 0; n < dimension; n++)
      {
        Sai1 << C[i][m][n] << " ";
      }
      Sai1 << std::endl;
    }
    Sai1 << " E= " << std::endl;
    for (int m = 0; m < dimension * dimension; m++)
    {
      Sai1 << E[i][m] << " ";
    }
    Sai1 << std::endl;
    Sai1 << " S= " << std::endl;
    for (int m = 0; m < dimension * dimension; m++)
    {
      Sai1 << S[i][m] << " ";
    }
    Sai1 << std::endl;
    Sai1 << " CC= " << std::endl;
    for (int m = 0; m < dimension * dimension; m++)
    {
      for (int n = 0; n < dimension * dimension; n++)
      {
        Sai1 << CC[i][m][n] << " ";
      }
      Sai1 << std::endl;
    }
    Sai1 << " fint= ";
    for (int j = 0; j < dimension * NNo; j++)
    {
      Sai1 << fint[i][j] << " ";
    }
    Sai1 << std::endl;
    Sai1 << " h= " << std::endl;
    for (int a = 0; a < dimension * NNo; a++)
    {
      for (int j = 0; j < dimension * NNo; j++)
      {
        Sai1 << h[i][a][j] << " ";
      }
      Sai1 << std::endl;
    }
    Sai1 << " ue[" << i << "]=" << ue[i] << " J0=" << J0[i] << " W[i]=" << r_Quad()->r_W(i) << std::endl;
    Sai1 << std::endl;
  }
  Sai1 << " Ue= " << Ue << std::endl;
  Sai1 << " Fint= " << std::endl;
  for (int i = 0; i < dimension * NNo; i++)
  {
    Sai1 << Fint[i] << " ";
  }
  Sai1 << std::endl;
  Sai1 << " H= " << std::endl;
  for (int i = 0; i < dimension * NNo; i++)
  {
    for (int j = 0; j < dimension * NNo; j++)
    {
      Sai1 << Hessiana[i][j] << " ";
    }
    Sai1 << std::endl;
  }
  Sai1 << " M= " << std::endl;
  for (int i = 0; i < dimension * NNo; i++)
  {
    for (int j = 0; j < dimension * NNo; j++)
    {
      Sai1 << Mass[i][j] << " ";
    }
    Sai1 << std::endl;
  }
  Sai1 << " MRP= " << std::endl;
  for (int i = 0; i < NNo; i++)
  {
    for (int j = 0; j < NNo; j++)
    {
      Sai1 << MRP[i][j] << " ";
    }
    Sai1 << std::endl;
  }

  Sai1 << std::endl;
  Sai1 << std::endl;
  Sai1.close();
}
