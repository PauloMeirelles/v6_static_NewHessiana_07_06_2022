//#pragma hdrstop
//#pragma argsused


//=============================
//     F5 -> Debbug
//     CTRL+G5 -> Compilar sem debbugar
//
//    obs: qualquer alteração no código, seguir as seguintes etapas:
//    1 - salvar os arquivos alterados;
//    2 - entrar na pasta build pelo terminal "cd build" e usar o comando "make" para refazer 
//    os arquivos .o e o executvel;
//    ARQUIVO VERSIONADO A SER UTILIZADO: /home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022
//    CAMINHO DE ENTRADA: config.txt referente ao caminho acima
//=================================

#include "main.h"

static char help[]="Solver Petsc";

int main(int argc, char *argv[])
{
    PetscInitialize(&argc, &argv, (char*)0, help);   //-> MODELO PARA MPI
//=====================================================
/*
    V2D A;
    A.resize(5,V1D(5,0.0));
    V1D b;
    b.resize(5,0.0);
    for(int i=0; i<5; i++)
    {
      for (int j=0; j<5; j++)
      {
        if(i!=j)
        A[i][j]=2;
        else
        A[i][j]=1;
      }
    }

    for(int i=0; i<5; i++)
    {
      b[i]=i+1;
    }

//====== ANYEIGEN ============
    V2D anyEigen;
    anyEigen.resize(5,V1D(5,0.0));
    anyEigen=invAnyEIGEN(A); std::cout << "AnyEigen" << "\n";
    for(int i=0; i<5; i++)
    {
      std::cout << "[";
      for (int j=0; j<5; j++)
      {
          std::cout << anyEigen[i][j] << ", ";
      }
            std::cout << "]" << "\n";
    } 

//====== invEigen ============
    V2D invAA;
    invAA.resize(5,V1D(5,0.0));
    invAA=invEIGEN(A); std::cout << "invEigen" << "\n";
    for(int i=0; i<5; i++)
    {
      std::cout << "[";
      for (int j=0; j<5; j++)
      {
          std::cout << invAA[i][j] << ", ";
      }
            std::cout << "]" << "\n";
    } 

//===== EINGEN SOLVER ==========
    V1D xES;
    xES.resize(5,0.0);
    xES=EIGENsolver(argc,argv,A,b); std::cout << "Eigen Solver" << "\n";
    for(int i=0; i<5; i++)
    {
      std::cout << "[";
      std::cout << xES[i] << ", ";
      std::cout << "]" << "\n";
    } 

//===== KSP SOLVER ==========
    V1D xKSP;
    xKSP.resize(5,0.0);
    xKSP=KSPresolve(argc,argv,A,b); std::cout << "KSP Solver" << "\n";
    for(int i=0; i<5; i++)
    {
      std::cout << "[";
      std::cout << xKSP[i] << ", ";
      std::cout << "]" << "\n";
    } 

//===== PETSC SOLVER ==========
    V1D xPetsc;
    Mat pA;
    MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                     5, 5,100,NULL,300,NULL,&pA);
    xPetsc.resize(5,0.0);
    for (int i = 0; i < 5; i++)
    {
      for (int j = 0; j < 5; j++)
      {
      const PetscInt* p_i=& i;
      const PetscInt* p_j=& j;
      const PetscScalar* p_A=&A[i][j]; //std::cout << "Hpetsc [" << *p_He << "]\n";
      MatSetValues(pA, 1, p_i, 1, p_j, p_A, ADD_VALUES);  //std::cout << H_static(p_i)(p_j) << "]\n";
      }
    }

    MatAssemblyBegin(pA, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(pA, MAT_FINAL_ASSEMBLY);
    MatView(pA, PETSC_VIEWER_STDOUT_WORLD); 
    xPetsc=PETScSolver(argc,argv,pA,b); std::cout << "Petsc Solver" << "\n";
    for(int i=0; i<5; i++)
    {
      std::cout << "[";
      std::cout << xPetsc[i] << ", ";
      std::cout << "]" << "\n";
    }
    */

  //BREAK

    clock_t Ticks[2];
    Ticks[0] = clock();

    int k = 0;
    std::string EndCON = "/home/oem/PRG/CPP/Paulo/MEFPoro/v6_static_NewHessiana_07_06_2022/config.txt";
    std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(7);
    tvNo No;
    tvEl El;
    tvMa Ma;
    tvQuadratura Qua;
    tvTempo vTempo;
    tvAPPMETHOD vAppMethod;
    CCo Co(EndCON);
    int ngl;
    // Apaga os arquivos a serem preenchidos com os relatorios...
    std::string NArq1;
    NArq1.append(Co.r_ArqREL()); // NArq1.append("resultados/relatorio.txt");
    std::ofstream Sai1;
    Sai1.open(NArq1.c_str()); ///< numero de nos e matriz a ser invertida
    Sai1.close();
    Le_Tempo(vTempo,Co.r_ArqTEMPO());
    Le_AppMethod(vAppMethod,Co.r_ArqMETHOD());
    Le_Quadraturas(Qua, Co.r_ArqQUA()); // ALTERAR O ENDEREÇO O CONTIDO NO PRIMEIRO OBJ DA CLASSE NO VETOR CRIADO <QUADRATURA E NUMERO DA QUADRATURA> INDIRETAMENTE COM UM PUSHBACK
    std::cout << Co.r_ArqNOS() << std::endl;
    Le_Nos(No, Co.r_ArqNOS());
    Le_CC(No, Co.r_ArqCC(), ngl);
    Le_P(No, Co.r_ArqP());
    Ma = Le_Mat(Co.r_ArqMAT());
    Le_Elementos(El, No, Ma, Qua, Co);

    //========== CONSTRUCT MINIMOS QUADRADOS GLOBAIS ====================
    int NNQ = 0;
    int tm = El.size();
    std::cout << "Tamanho do vetor de elementos: " << tm << "\n";
    for (int i = 0; i < El.size(); i++)
    {
        NNQ += El[i].r_Quad()->r_NumberQuad();
    }
    int NNP = 0; //!!!!!!!cuidado esta variavel e tudo que ela implica tem que ser seriamente repensada
    /// especialmente para os casos em que existem n�s declarados que nao sao usados na estrutura....
    for (int i = 0; i < No.size(); i++)
    {
      NNP += No[i].r_FP(); // NOS QUE FAZEM PARTE
    }
    // std::cout << "NNQ:" << NNQ << " && NNP:" << NNP << std::endl;
    int nj = 0;
    // ==========================================================================================

   int ApproxMethod=vAppMethod[0].r_Method();
   int OutputFile=vAppMethod[0].r_OutputFile();

V2D iMPHI;
if (ApproxMethod==0) 
{
    //======== Metodo Minimos Quadrados Globais (PseudoInversa) =========================
    V2D MPHI;
    MPHI.resize(NNQ, V1D(NNP, 0.0));
    for (int i = 0; i < El.size(); i++)
    {
        for (int j = 0; j < El[i].r_Quad()->r_NumberQuad(); j++)
        {
            for (int k = 0; k < El[i].r_NNo(); k++)
            {
                MPHI[nj][El[i].r_No(k)->r_NumberNo()] = El[i].r_PHI(j, k);
            }
            nj = nj + 1;
        }
    }
    int n = MPHI.size();
    int m = MPHI[1].size();
    //std::cout << "VALOR DE N: " << n << "\n";
    //std::cout << "VALOR DE M: " << m << "\n";


    // iMPHI.resize(MPHI[0].size(),V1D (MPHI.size(),0.0));
    /*
    std::cout << "MPHI"
              << "\n";
    for (int t = 0; t < n; t++)
    {
        std::cout << "[";
        for (int q = 0; q < m; q++)
        {
            std::cout << MPHI[t][q] << ",";
        }
        std::cout << "]"
                  << "\n";
    }
    std::cout << "FIM!"
              << "\n";
    */

    iMPHI = pseudoInverse(MPHI);
}
//=============================================================

    std::vector<std::vector<double>> Fint;
    Fint.resize(ngl, std::vector<double>(1, 0.0));
    std::vector<std::vector<double>> Fext;
    Fext.resize(ngl, std::vector<double>(1, 0.0));
    std::vector<double> _g;
    _g.resize(ngl, 0.0);
    std::vector<double> Dy;
    Dy.resize(ngl, 0.0);
    std::vector<std::vector<double>> H;
    H.resize(ngl, std::vector<double>(ngl, 0.0));
    std::vector<std::vector<double>> M;
    M.resize(ngl, std::vector<double>(ngl, 0.0));

    // ============== PARTE DO TESTE =======
    Mat  H_static, M_global;
    // std::vector<double> Dy_;
    // Dy_.resize(ngl, 0.0);
    

    //MatCreate(PETSC_COMM_WORLD, &H_static);
    //MatSetSizes(H_static, PETSC_DECIDE, PETSC_DECIDE, ngl, ngl);
    //MatSetFromOptions(H_static);
    //MatSetUp(H_static);

    //MatCreate(PETSC_COMM_WORLD, &M_global);
    //MatSetSizes(M_global, PETSC_DECIDE, PETSC_DECIDE, ngl, ngl);
    //MatSetFromOptions(M_global);
    //MatSetUp(M_global);
    // ===============================

    double Nx;
    double NDy;
    //double NDy_teste;
    clock_t TicksSolver[2];
    
    Nx = 0;
    for (int i = 0; i < No.size(); i++)
    {
        for (int j = 0; j < No[i].r_dimension(); j++) //Aqui esta lendo nas 3 direções (ALTERAR PARA Dimension)
        {
            if (No[i].r_CondCountour(j) == 0)
            {
                Nx += pow(No[i].r_X(0, j), 2.0);
            }
        }
    }
    Nx = pow(Nx, 0.5);
    double Y=0;
    int ite;
    ite = 0;

    
  for(int fs=0; fs<vTempo[0].r_ForceStep(); fs++)
  {
    double tol=vTempo[0].r_Tolerance();
    double FStep=vTempo[0].r_ForceStep();
    double nFS=1;
    nFS+=fs;
    std::cout << "===========================================================" << "\n";
    std::cout << "PASSO DE CARGA ATUAL:" << nFS << std::endl;
    std::cout << "===========================================================" << "\n";
    do
    { 
      /*
      for (int i = 0; i < El.size(); i++)
      {
        El[i].ite();
        for (int a = 0; a < El[i].r_NNo(); a++)
        {
          for (int z = 0; z < 2; z++)
          {
            if (El[i].r_No(a)->r_Adress(z) >= 0)
            {
                Fint[El[i].r_No(a)->r_Adress(z)][0] += El[i].r_Fint(El[i].r_dimension() * a + z);
                Fext[El[i].r_No(a)->r_Adress(z)][0] = ((El[i].r_No(a)->r_PConcentrated(z))*(nFS/FStep));
                _g[El[i].r_No(a)->r_Adress(z)] = Fext[El[i].r_No(a)->r_Adress(z)][0] - Fint[El[i].r_No(a)->r_Adress(z)][0];
                for (int l = 0; l < El[i].r_NNo(); l++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        if (El[i].r_No(l)->r_Adress(k) >= 0)
                        {
                          H[El[i].r_No(a)->r_Adress(z)][El[i].r_No(l)->r_Adress(k)] += El[i].r_Hessiana(El[i].r_dimension() * a + z, El[i].r_dimension() * l + k);
                          M[El[i].r_No(a)->r_Adress(z)][El[i].r_No(l)->r_Adress(k)] += El[i].r_Mass(El[i].r_dimension() * a + z, El[i].r_dimension() * l + k);
                        }
                    }
                }
            }
          }
        }
      } */

      // ====================== TESTE EM VARIAVEIS DO PETSC ===========
      //V1D g_;
      //g_=StaticEq(argc, argv, El, ngl, nFS, FStep);
      
     MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                     ngl, ngl,100,NULL,300,NULL,&H_static);
     MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                     ngl, ngl,100,NULL,300,NULL,&M_global);
      // Construindo aqui fora por enquanto
      for (int i = 0; i < El.size(); i++)
      {
        std::cout << "Elemento: " << i << "\n";
        El[i].ite();
        for (int a = 0; a < El[i].r_NNo(); a++)
        {
          for (int z = 0; z < 2; z++)
          {
            if (El[i].r_No(a)->r_Adress(z) >= 0)
            {
                Fint[El[i].r_No(a)->r_Adress(z)][0] += El[i].r_Fint(El[i].r_dimension() * a + z);
                Fext[El[i].r_No(a)->r_Adress(z)][0] = ((El[i].r_No(a)->r_PConcentrated(z))*(nFS/FStep));
                _g[El[i].r_No(a)->r_Adress(z)] = (Fint[El[i].r_No(a)->r_Adress(z)][0]-Fext[El[i].r_No(a)->r_Adress(z)][0])*-1;
                for (int l = 0; l < El[i].r_NNo(); l++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        if (El[i].r_No(l)->r_Adress(k) >= 0)
                        {
                          int ii=0; ii=El[i].r_No(a)->r_Adress(z);
                          int jj=0; jj=El[i].r_No(l)->r_Adress(k);
                          int c=0; c=El[i].r_dimension() * a + z; //std::cout << c <<"\n";
                          int b=0; b=El[i].r_dimension() * l + k; //std::cout << c <<"\n";
                          std::cout << "Hlocal [" << c << "]" << "[" << b << "] " << El[i].r_Hessiana(El[i].r_dimension() * a + z, El[i].r_dimension() * l + k)<< "\n";
                          std::cout << "Hglobal [" << El[i].r_No(a)->r_Adress(z) << "]" << "[" << El[i].r_No(l)->r_Adress(k) << "] " << "\n";
                          H[El[i].r_No(a)->r_Adress(z)][El[i].r_No(l)->r_Adress(k)] += El[i].r_Hessiana(El[i].r_dimension() * a + z, El[i].r_dimension() * l + k);
                          M[El[i].r_No(a)->r_Adress(z)][El[i].r_No(l)->r_Adress(k)] += El[i].r_Mass(El[i].r_dimension() * a + z, El[i].r_dimension() * l + k);
                                      /*  for(int linha=0; linha <24; linha++)
                                       {
                                         std::cout << "Hessiana: Row "<< linha <<" [";
                                        for (int coluna = 0; coluna< 24; coluna++)
                                        {
                                        std::cout << H[linha][coluna]<< ",";
                                        }
                                         std::cout <<" ]" << "\n";

                                       }*/
                          //double H_sum=0.0; H_sum += El[i].r_Hessiana(El[i].r_dimension() * a + z, El[i].r_dimension() * l + k);
                          //double M_sum=0.0; M_sum += El[i].r_Mass(El[i].r_dimension() * a + z, El[i].r_dimension() * l + k);
                          //std::cout << "H [" << ii << "]" << "[" << jj << "] " << H[El[i].r_No(a)->r_Adress(z)][El[i].r_No(l)->r_Adress(k)] << "} ";
                          const PetscInt* p_i=& ii;
                          const PetscInt* p_j=& jj;
                          const PetscScalar* p_He=&H[ii][jj]; //std::cout << "Hpetsc [" << *p_He << "]\n";
                          const PetscScalar* p_Ma=&H[ii][jj];
                          MatSetValues(H_static, 1, p_i, 1, p_j, p_He, INSERT_VALUES);  //std::cout << H_static(p_i)(p_j) << "]\n";
                          MatSetValues(M_global, 1, p_i, 1, p_i, p_Ma , INSERT_VALUES);
                        }
                    }
                }
            }
          }
        }
      }
      // ====================== FIM TESTE EM VARIAVEIS DO PETSC =======
      MatAssemblyBegin(H_static, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(H_static, MAT_FINAL_ASSEMBLY);
      MatView(H_static, PETSC_VIEWER_STDOUT_WORLD); 




        // ====== VISUALIZAÇAO DE RESULTADOS
           /*
           std::cout << "Matriz H Iteracao=" << ite << std::endl;
           std::cout << "MatrixForm[HH={";
           for (int i=0; i<H.size(); i++)
           {
               std::cout << "{";
               for (int k=0; k<H[0].size(); k++)
               {
                   std::cout << H[i][k];
                   if (k<H[0].size()-1) std::cout << ",";
               }
               std::cout << "}";
               if (i<H.size()-1) std::cout << ",";
           }
           std::cout << std::endl;
        */
        /* // ===== AQUI PODE ABRIR PARA VISUALIZAR
       std::cout << "Matriz M Iteracao=" << ite << std::endl;
       std::cout << "MatrixForm[MM={";
       for (int i=0; i<M.size(); i++)
       {
       std::cout << "{";
       for (int k=0; k<M[0].size(); k++)
       {
       std::cout << M[i][k];
       if (k<M[0].size()-1) std::cout << ",";
       }
       std::cout << "}";
       if (i<M.size()-1) std::cout << ",";
       }
        std::cout << "}]";
        std::cout << std::endl;
        */

        /*
        std::cout << "Vetor F Iteracao=" << ite << std::endl;
        std::cout << "MatrixForm[FF={";
        for (int i=0; i<_g.size(); i++)
        {
            std::cout << _g[i];
            if (i<_g.size()-1) std::cout << ",";
        }
        std::cout << "}]";
        std::cout << std::endl;
         */

        Sai1.open(NArq1.c_str(), std::fstream::in | std::fstream::out | std::fstream::app); ///< Dados Calculados

        //==================== RESOLUÇÃO DO SISTEMA LINEAR ==================================
        //TicksSolver[0]=0;TicksSolver[1]=0;
        TicksSolver[0] = clock();
        //Dy=PETScSolver(argc, argv, H_static, _g);
        //Dy=KSPresolve(argc, argv, H, _g);
        Dy=EIGENsolver(argc, argv, H, _g);
        TicksSolver[1] = clock();
        double TimeSolver = (TicksSolver[1] - TicksSolver[0]) *1000.0 / CLOCKS_PER_SEC;
        std::cout << "\tIteraçao " << ite << " - Tempo gasto para construção do sistema e Resoluçao :" << std::setprecision(7) << TimeSolver << "ms." << "\n";
        //Dy = MUMPSresolve(argc, argv, H, _g); //--> NAO OK!
        //Dy= MUMPSsolver(argc, argv, H, _g);


        /*
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "============================================" << std::endl;
        std::cout << "Valores de Dy:"
                  << "\n";
        for (int i = 0; i < ngl; i++)
        {
            std::cout << Dy[i] << " ";
        }
        std::cout << std::endl;
        std::cout << std::endl;
        */
        
        NDy = 0;
        for (int i = 0; i < Dy.size(); i++)
        {
          NDy += pow(Dy[i], 2.0);
        }
        NDy = pow(NDy, 0.5);

        // NDy_teste = 0;
        // for (int i = 0; i < Dy.size(); i++)
        // {
        //   NDy_teste += fabs(Dy[i]);
        // }

        Sai1 << " Dy = ";
        for (int i = 0; i < Dy.size(); i++)
        {
          Sai1 << Dy[i] << " ";
        }
        Sai1 << std::endl;
        Sai1.close();
        // ========== ANTERIOR ESTATICO FORA DO MAT PETSC ==========
        /*
        Sai1 << " Dy = ";
        for (int i = 0; i < Dy.size(); i++)
        {
          Sai1 << Dy[i] << " ";
        }
        Sai1 << std::endl;
        Sai1.close();
        */
        // ========== ANTERIOR ESTATICO FORA DO MAT PETSC ==========
        for (int i = 0; i < ngl; i++)
        {
          Fint[i][0] = 0.0;
            for (int j = 0; j < ngl; j++)
            {
              H[i][j] = 0.0;
              M[i][j] = 0.0;
            }
        }
        k = 0;
        Y=0;
        for (int i = 0; i < No.size(); i++)
        {
            for (int j = 0; j < 2; j++)
            {
                if (No[i].r_CondCountour(j) == 0)
                {
                  k = k + 1; ///++K atualiza antes k++ atualiza depois;
                  Y = No[i].r_X(1, j) + Dy[k - 1];
                  No[i].w_X(1, j, Y);
                }
            }
            // std::cout << i << " " << No[i].r_X(1,0) << " " << No[i].r_X(1,1) << " " << No[i].r_X(1,2) << std::endl;
        }
        for(int u=0; u< ngl; u++)
        {
        Dy[u]=0.0;
        }
      if(NDy / Nx < tol) 
      {
        // Metodo Minimos Quadrados Globais =========================
        if (ApproxMethod==0) 
        {
            std::cout << "Condição NDy/Nx válida - Calculo com Minimos Quadrados (PseudoInversa)" << std::endl;

            for (int k1 = 0; k1 < 2; k1++)
            {
                for (int k2 = 0; k2 < 2; k2++)
                {
                    V1D GRA;
                    GRA.resize(NNQ, 0.0);
                    V1D GRB;
                    GRB.resize(NNP, 0.0);
                    int a = 0;
                    for (int i = 0; i < El.size(); i++)
                    {
                        for (int j = 0; j < El[i].r_Quad()->r_NumberQuad(); j++)
                        {
                          GRA[a] = El[i].r_E(j, 2 * k1 + k2);
                          a = a + 1;
                        }
                    }
                    GRB = MV(iMPHI, GRA);
                    for (int i = 0; i < No.size(); i++)
                    {
                       No[i].w_no_E(2 * k1 + k2, GRB[i]);
                    }
                    for (int i = 0; i < El.size(); i++)
                    {
                        for (int j = 0; j < El[i].r_NNo(); j++)
                        {
                          El[i].w_no_E(j, 2 * k1 + k2, El[i].r_No(j)->r_no_E(2 * k1 + k2));
                        }
                    }
                }
            }
            for (int k1 = 0; k1 < 2; k1++)
            {
                for (int k2 = 0; k2 < 2; k2++)
                {
                    V1D GRA;
                    GRA.resize(NNQ, 0.0);
                    V1D GRB;
                    GRB.resize(NNP, 0.0);
                    int a = 0;
                    for (int i = 0; i < El.size(); i++)
                    {
                        for (int j = 0; j < El[i].r_Quad()->r_NumberQuad(); j++)
                        {
                          GRA[a] = El[i].r_S(j, 2 * k1 + k2);
                          a = a + 1;
                        }
                    }
                    GRB = MV(iMPHI, GRA);
                    for (int i = 0; i < No.size(); i++)
                    {
                      No[i].w_no_S(2 * k1 + k2, GRB[i]);
                    }
                    for (int i = 0; i < El.size(); i++)
                    {
                        for (int j = 0; j < El[i].r_NNo(); j++)
                        {
                          El[i].w_no_S(j, 2 * k1 + k2, El[i].r_No(j)->r_no_S(2 * k1 + k2));
                        }
                    }
                }
            }
            V1D GRAa;
            GRAa.resize(NNQ, 0.0);
            V1D GRBa;
            GRBa.resize(NNP, 0.0);
            int aa = 0;
            for (int i = 0; i < El.size(); i++)
            {
                for (int j = 0; j < El[i].r_Quad()->r_NumberQuad(); j++)
                {
                  GRAa[aa] = El[i].r_ue(j);
                  aa = aa + 1;
                }
            }
            GRBa = MV(iMPHI, GRAa);
            for (int i = 0; i < No.size(); i++)
            {
              No[i].w_no_ue(GRBa[i]);
            }
            for (int i = 0; i < El.size(); i++)
            {
                for (int j = 0; j < El[i].r_NNo(); j++)
                {
                  El[i].w_no_ue(j, El[i].r_No(j)->r_no_ue());
                }
            }
        }
        else // ========================= COLOCAR CALC SIGMA ======== (ACIMA)
        {
          std::cout << "Condição NDy/Nx válida - Calculo com Nó a Nó + Sem Media no contorno" << std::endl;
          for (int i=0; i<No.size(); i++)
          {
            for (int k=0; k<4; k++) 
            {
              No[i].w_no_E(k,0.0);
              No[i].w_no_S(k,0.0);
              No[i].w_no_Sig(k,0.0);
                for (int j=0; j<No[i].sizeElm(); j++)
                {
                  No[i].w_no_E (k,No[i].r_no_E (k)+El[No[i].r_Elm(j)].r_no_E (No[i].r_no_Elm(j))[k]);
                  No[i].w_no_S (k,No[i].r_no_S (k)+El[No[i].r_Elm(j)].r_no_S (No[i].r_no_Elm(j))[k]);
                  No[i].w_no_Sig (k,No[i].r_no_Sig (k)+El[No[i].r_Elm(j)].r_no_Sig (No[i].r_no_Elm(j))[k]);
                }
              No[i].w_no_E (k,No[i].r_no_E (k)/No[i].sizeElm());
              No[i].w_no_S (k,No[i].r_no_S (k)/No[i].sizeElm());
              No[i].w_no_Sig (k,No[i].r_no_Sig (k)/No[i].sizeElm());
            }
            No[i].w_no_ue(0.0);    
            for (int j=0; j<No[i].sizeElm(); j++)
            {
              No[i].w_no_ue(No[i].r_no_ue()+El[No[i].r_Elm(j)].r_no_ue(No[i].r_no_Elm(j)));
            }
            No[i].w_no_ue(No[i].r_no_ue()/No[i].sizeElm());
          }
        }
        //====================== FIM DO METODO DE MINIMOS QUADRADOS GLOBAIS =====================================
        
        // ==== GERAÇãO DO ARQUIVO DE SAIDA ====
        std::cout << "Iteração " << ite << " convergiu - Gerando arquivo de saida!" << std::endl;
        {
        if(OutputFile==0)
        {
          Gera_VTK_69(Co, No, El, fs);
        }
          else if(OutputFile==1)
          {
            Gera_VTK_69_NI(Co, No, El, fs);     
          }
            else
            {
                Gera_VTU_69(Co, No, El, fs);
            }
        }
        //
        //Gera_VTU_69(Co, No, El, DeltaP);
      } 
        //std::cout << "Iteração número:" << ite << std::endl;
        // std::cout << ite << std::endl;
        ite += 1;

    } while (NDy / Nx > tol); //FIM DO FULL-NEWTON-RAPHSON COM CONTROLE DE POSIÇãO
  } //FIM DO INCREMENTO DE CARGA
    
    std::cout << "\n";
    std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
              << "\n";
    std::cout << "Processo Finalizado ! PETSc Finalize invocado." << std::endl;
    std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
              << "\n";
    PetscFinalize();
    Ticks[1] = clock();
    double Tempo = (Ticks[1] - Ticks[0]) *1000.0 / CLOCKS_PER_SEC;
    std::cout << "Tempo gasto:" << std::setprecision(7) << Tempo << "ms.";
    getchar(); 
    return 0;
}

void Gera_VTU_69(CCo &Co, tvNo &No, tvEl &El, const int &ite)
{
    std::cout << "Gerando arquivos de saída formato vtu (Interpolado) .."  << "\n";
    std::vector<std::vector<int>> ind;
    ind.resize(2, std::vector<int>(2, 0.0));
    ind[0][0] = 0;
    ind[0][1] = 2;
    ind[1][0] = 3;
    ind[1][1] = 1;
    std::string NArq2;
    NArq2.append(Co.r_ArqVTU());
    NArq2.insert(NArq2.find(".vtu"), intToString(ite));
    std::ofstream Sai2;
    Sai2.open(NArq2.c_str());
    Sai2.close();
    Sai2.open(NArq2.c_str(), std::fstream::in | std::fstream::out | std::fstream::app); ///< Dados Calculado
    // std::stringstream text;
	// text << "output" << ite << ".vtu";
	// std::ofstream Sai2(text.str());
    

    //Header
	Sai2 << "<?xml version=\"1.0\"?>" << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">" << "\n"
		 << "  <UnstructuredGrid>" << "\n"
         << "  <Piece NumberOfPoints=\"" << No.size()
         << "\"  NumberOfCells=\"" << El.size()
         << "\">" << "\n";
	//Nodal Coordinates
	Sai2 << "    <Points>" << "\n"
         << "      <DataArray type=\"Float64\" "
         << "NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_X(1, 0) << " " << No[i].r_X(1, 1) << " " << 0.0 << std::endl;
    }
    //Sai2 << std::endl;
    Sai2 << "      </DataArray>" << "\n"
         << "    </Points>" << "\n";
	//Element Connectivity
	Sai2 << "    <Cells>" << "\n"
         << "      <DataArray type=\"Int32\" "
         << "Name=\"connectivity\" format=\"ascii\">" << "\n";
    for (int i = 0; i < El.size(); i++)
    {
        //Sai2 << "10";
        for (int j = 0; j < El[i].r_NNo(); j++)
        {
            Sai2 << " " << El[i].r_No(j)->r_NumberNo();
        }
        Sai2 << std::endl;
    }

    Sai2 << "      </DataArray>" << "\n";
	//offsets
	Sai2 << "      <DataArray type=\"Int32\""
		 << " Name=\"offsets\" format=\"ascii\">" << "\n";
    int aux=0;
	for (int i=0; i< El.size(); i++)
	{
		int n = El[i].r_NNo();
		aux += n;
		Sai2 << aux << "\n";
	}
	Sai2 << "      </DataArray>" << "\n";
	//elements type
	Sai2 << "      <DataArray type=\"UInt8\" Name=\"types\" "
		 << "format=\"ascii\">" << "\n";
    for (int i = 0; i < El.size(); i++)
    {
        Sai2 << "69" << std::endl;
    }
    //Sai2 << std::endl;
	Sai2 << "      </DataArray>" << "\n"
		 << "    </Cells>" << "\n";
	//nodal results Displacement
	Sai2 << "    <PointData>" <<"\n";
	Sai2 << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
		<< "Name=\"Displacement\" format=\"ascii\">" << "\n";
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_X(1, 0) - No[i].r_X(0, 0) << " " << No[i].r_X(1, 1) - No[i].r_X(0, 1) << std::endl;
    }
    //Sai2 << std::endl;

   
    //nodal results Tensor S
	Sai2 << "      </DataArray>" << "\n";
	Sai2 << "      <DataArray type=\"Float64\" NumberOfComponents=\"9\" "
		<< "Name=\"Tensor S\" format=\"ascii\">" << "\n";
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_no_S(ind[0][0]) << " " << No[i].r_no_S(ind[0][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << No[i].r_no_S(ind[1][0]) << " " << No[i].r_no_S(ind[1][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << "0.00"
             << " "
             << "0.00"
             << " "
             << " 0.00" << std::endl;
        Sai2 << std::endl;
    }
    //Sai2 << std::endl;

    //nodal results Tensor E
	Sai2 << "      </DataArray>" << "\n";
	Sai2 << "      <DataArray type=\"Float64\" NumberOfComponents=\"9\" "
		<< "Name=\"Tensor E\" format=\"ascii\">" << "\n";
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_no_E(ind[0][0]) << " " << No[i].r_no_E(ind[0][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << No[i].r_no_E(ind[1][0]) << " " << No[i].r_no_E(ind[1][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << "0.00"
             << " "
             << "0.00"
             << " "
             << " 0.00" << std::endl;
        Sai2 << std::endl;
    }
      
    //Sai2 << std::endl;

    //nodal results Tensor Sigma
	Sai2 << "      </DataArray>" << "\n";
	Sai2 << "      <DataArray type=\"Float64\" NumberOfComponents=\"9\" "
		<< "Name=\"Sigma\" format=\"ascii\">" << "\n";
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_no_Sig(ind[0][0]) << " " << No[i].r_no_Sig(ind[0][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << No[i].r_no_Sig(ind[1][0]) << " " << No[i].r_no_Sig(ind[1][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << "0.00"
             << " "
             << "0.00"
             << " "
             << " 0.00" << std::endl;
        Sai2 << std::endl;
    }
    //Sai2 << std::endl;

    //nodal results Tensor Sigma
	Sai2 << "      </DataArray>" << "\n";
	Sai2 << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
		<< "Name=\"ue\" format=\"ascii\">" << "\n";
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_no_ue() << std::endl;
    }
    //Sai2 << std::endl;
        
	// Sai2 << "      </DataArray> " << "\n";
	// Sai2 << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
	// 	<< "Name=\"Velocity\" format=\"ascii\">" << "\n";

	// for (Node* n: nodes_)
	// {
	// 	file << n->getCurrentVelocity()(0) << " " << n->getCurrentVelocity()(1) << "\n";
	// }
	// file << "      </DataArray> " << "\n";
	// file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
	// 	<< "Name=\"Acceleration\" format=\"ascii\">" << "\n";

	// for (Node* n: nodes_)
	// {
	// 	file << n->getCurrentAcceleration()(0) << " " << n->getCurrentAcceleration()(1) << "\n";
	// }
    
	Sai2 << "      </DataArray> " << "\n";
	Sai2 << "    </PointData>" << "\n";
	//Elemental results
	Sai2 << "    <CellData>" << "\n";

	Sai2 << "    </CellData>" << "\n";
	//Footnote
	Sai2 << "  </Piece>" << "\n"
		<< "  </UnstructuredGrid>" << "\n"
		<< "</VTKFile>" << "\n";
    Sai2.close();
}

void Gera_VTK_69(CCo &Co, tvNo &No, tvEl &El, const int &ite)
{
    std::cout << "Gerando arquivos de saída formato vtk (Interpolado) .."
              << "\n";
    std::vector<std::vector<int>> ind;
    ind.resize(2, std::vector<int>(2, 0.0));
    ind[0][0] = 0;
    ind[0][1] = 2;
    ind[1][0] = 3;
    ind[1][1] = 1;
    std::string NArq2;
    NArq2.append(Co.r_ArqVTK());
    NArq2.insert(NArq2.find(".vtk"), intToString(ite));
    std::ofstream Sai2;
    Sai2.open(NArq2.c_str());
    Sai2.close();
    Sai2.open(NArq2.c_str(), std::fstream::in | std::fstream::out | std::fstream::app); ///< Dados Calculado
    Sai2 << "# vtk DataFile Version 2.0" << std::endl;
    Sai2 << "FEM, Created by Wutzow" << std::endl;
    Sai2 << "ASCII" << std::endl;
    Sai2 << "DATASET UNSTRUCTURED_GRID" << std::endl;
    Sai2 << "POINTS " << No.size() << " Float64" << std::endl;
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_X(1, 0) << " " << No[i].r_X(1, 1) << " " << 0.0 << std::endl;
    }
    Sai2 << std::endl;
    Sai2 << "CELLS " << El.size() << " " << El.size() * 11 << std::endl;
    for (int i = 0; i < El.size(); i++)
    {
        Sai2 << "10";
        for (int j = 0; j < El[i].r_NNo(); j++)
        {
            Sai2 << " " << El[i].r_No(j)->r_NumberNo();
        }
        Sai2 << std::endl;
    }
    Sai2 << std::endl;
    Sai2 << "CELL_TYPES " << El.size() << std::endl;
    for (int i = 0; i < El.size(); i++)
    {
        Sai2 << "69" << std::endl;
    }
    Sai2 << std::endl;
    Sai2 << "CELL_DATA " << El.size() << std::endl;
    Sai2 << "POINT_DATA " << No.size() << std::endl;
    Sai2 << "VECTORS Y"
         << " Float64" << std::endl;
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_X(1, 0) << " " << No[i].r_X(1, 1) << " 0.0" << std::endl;
    }
    Sai2 << std::endl;
    Sai2 << "CELL_DATA " << El.size() << std::endl;
    Sai2 << "POINT_DATA " << No.size() << std::endl;
    Sai2 << "VECTORS U"
         << " Float64" << std::endl;
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_X(1, 0) - No[i].r_X(0, 0) << " " << No[i].r_X(1, 1) - No[i].r_X(0, 1) << " 0.0" << std::endl;
    }
    Sai2 << std::endl;

    Sai2 << "CELL_DATA " << El.size() << std::endl;
    Sai2 << "POINT_DATA " << No.size() << std::endl;
    Sai2 << "TENSORS E Float64" << std::endl;
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_no_E(ind[0][0]) << " " << No[i].r_no_E(ind[0][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << No[i].r_no_E(ind[1][0]) << " " << No[i].r_no_E(ind[1][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << "0.00"
             << " "
             << "0.00"
             << " "
             << " 0.00" << std::endl;
        Sai2 << std::endl;
    }
    Sai2 << std::endl;

    Sai2 << "CELL_DATA " << El.size() << std::endl;
    Sai2 << "POINT_DATA " << No.size() << std::endl;
    Sai2 << "TENSORS S Float64" << std::endl;
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_no_S(ind[0][0]) << " " << No[i].r_no_S(ind[0][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << No[i].r_no_S(ind[1][0]) << " " << No[i].r_no_S(ind[1][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << "0.00"
             << " "
             << "0.00"
             << " "
             << " 0.00" << std::endl;
        Sai2 << std::endl;
    }
    Sai2 << std::endl;

    Sai2 << "CELL_DATA " << El.size() << std::endl;
    Sai2 << "POINT_DATA " << No.size() << std::endl;
    Sai2 << "TENSORS Sigma Float64" << std::endl;
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_no_Sig(ind[0][0]) << " " << No[i].r_no_Sig(ind[0][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << No[i].r_no_Sig(ind[1][0]) << " " << No[i].r_no_Sig(ind[1][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << "0.00"
             << " "
             << "0.00"
             << " "
             << " 0.00" << std::endl;
        Sai2 << std::endl;
    }
    Sai2 << std::endl;

    Sai2 << "CELL_DATA " << El.size() << std::endl;
    Sai2 << "POINT_DATA " << No.size() << std::endl;
    Sai2 << "SCALARS ue float" << std::endl;
    Sai2 << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_no_ue() << std::endl;
    }
    Sai2 << std::endl;
    Sai2.close();
}



void Gera_VTK_69_NI(CCo &Co, tvNo &No, tvEl &El, const int &ite)
{
    std::cout << "Gerando arquivos de saída formato vtk (Nao interpolado) .."
              << "\n";
    std::vector<std::vector<int>> ind;
    ind.resize(2, std::vector<int>(2, 0.0));
    ind[0][0] = 0;
    ind[0][1] = 2;
    ind[1][0] = 3;
    ind[1][1] = 1;
    std::string NArq2;
    NArq2.append(Co.r_ArqVTK());
    NArq2.insert(NArq2.find(".vtk"), intToString(ite));
    std::ofstream Sai2;
    Sai2.open(NArq2.c_str());
    Sai2.close();
    Sai2.open(NArq2.c_str(), std::fstream::in | std::fstream::out | std::fstream::app); ///< Dados Calculado
    Sai2 << "# vtk DataFile Version 2.0" << std::endl;
    Sai2 << "FEM, Created by Wutzow" << std::endl;
    Sai2 << "ASCII" << std::endl;
    Sai2 << "DATASET UNSTRUCTURED_GRID" << std::endl;
    Sai2 << "POINTS " << 10*El.size() << " Float64" << std::endl;
    for (int i = 0; i < El.size(); i++)
    {
      for (int j = 0; j < 10; j++)
      { 
        Sai2 << El[i].r_No(j)->r_X(1, 0) << " " << El[i].r_No(j)->r_X(1, 1) << " " << 0.0 << std::endl;
      }
    }
    Sai2 << std::endl;
    Sai2 << "CELLS " << El.size() << " " << El.size() * 11 << std::endl;
    for (int i = 0; i < El.size(); i++)
    {
        Sai2 << "10";
        for (int j = 0; j < El[i].r_NNo(); j++)
        {
            Sai2 << " " << 10*i+j;
        }
        Sai2 << std::endl;
    }
    Sai2 << std::endl;
    Sai2 << "CELL_TYPES " << El.size() << std::endl;
    for (int i = 0; i < El.size(); i++)
    {
        Sai2 << "69" << std::endl;
    }
    Sai2 << std::endl;
    Sai2 << "CELL_DATA " << El.size() << std::endl;
    Sai2 << "POINT_DATA " << 10*El.size() << std::endl;
    Sai2 << "VECTORS Y"
         << " Float64" << std::endl;
    for (int i = 0; i < El.size(); i++)
    {
        for (int j = 0; j < El[i].r_NNo(); j++)
        {
          Sai2 << El[i].r_no_Y(j,0) << " " << El[i].r_no_Y(j,1) << " 0.0" << std::endl;
        }
    }
    Sai2 << std::endl;
    Sai2 << "CELL_DATA " << El.size() << std::endl;
    Sai2 << "POINT_DATA " << 10*El.size() << std::endl;
    Sai2 << "VECTORS U"
         << " Float64" << std::endl;
    for (int i = 0; i < El.size(); i++)
    {
        for (int j = 0; j < El[i].r_NNo(); j++)
        {
          Sai2 << El[i].r_no_U(j,0) << " " << El[i].r_no_U(j,1) << " 0.0" << std::endl;
        }
    }
    Sai2 << std::endl;

    Sai2 << "CELL_DATA " << El.size() << std::endl;
    Sai2 << "POINT_DATA " << 10*El.size() << std::endl;
    Sai2 << "TENSORS E Float64" << std::endl;
    for (int i = 0; i < El.size(); i++)
    {
        for (int j = 0; j < El[i].r_NNo(); j++)
        {
        Sai2 << El[i].r_no_E(j)[ind[0][0]] << " " << El[i].r_no_E(j)[ind[0][1]] << " "
             << " 0.00" << std::endl;
        Sai2 << El[i].r_no_E(j)[ind[1][0]] << " " << El[i].r_no_E(j)[ind[1][1]] << " "
             << " 0.00" << std::endl;
        Sai2 << "0.00"
             << " "
             << "0.00"
             << " "
             << " 0.00" << std::endl;
        Sai2 << std::endl;
        }
    }
    Sai2 << std::endl;

    Sai2 << "CELL_DATA " << El.size() << std::endl;
    Sai2 << "POINT_DATA " << No.size() << std::endl;
    Sai2 << "TENSORS Sigma Float64" << std::endl;
    for (int i = 0; i < No.size(); i++)
    {
        Sai2 << No[i].r_no_Sig(ind[0][0]) << " " << No[i].r_no_Sig(ind[0][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << No[i].r_no_Sig(ind[1][0]) << " " << No[i].r_no_Sig(ind[1][1]) << " "
             << " 0.00" << std::endl;
        Sai2 << "0.00"
             << " "
             << "0.00"
             << " "
             << " 0.00" << std::endl;
        Sai2 << std::endl;
    }
    Sai2 << std::endl;

    Sai2 << "CELL_DATA " << El.size() << std::endl;
    Sai2 << "POINT_DATA " << 10*El.size() << std::endl;
    Sai2 << "TENSORS S Float64" << std::endl;
    for (int i = 0; i < El.size(); i++)
    {
        for (int j = 0; j < El[i].r_NNo(); j++)
        {
        Sai2 << El[i].r_no_S(j)[ind[0][0]] << " " << El[i].r_no_S(j)[ind[0][1]] << " "
             << " 0.00" << std::endl;
        Sai2 << El[i].r_no_S(j)[ind[1][0]] << " " << El[i].r_no_S(j)[ind[1][1]] << " "
             << " 0.00" << std::endl;
        Sai2 << "0.00"
             << " "
             << "0.00"
             << " "
             << " 0.00" << std::endl;
        Sai2 << std::endl;
        }
    }
    Sai2 << std::endl;

    Sai2 << "CELL_DATA " << El.size() << std::endl;
    Sai2 << "POINT_DATA " << 10*El.size() << std::endl;
    Sai2 << "SCALARS ue float" << std::endl;
    Sai2 << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < El.size(); i++)
    {
        for (int j = 0; j < El[i].r_NNo(); j++)
        {
        Sai2 << El[i].r_no_ue(j) << std::endl;
        }
    }
    Sai2 << std::endl;
    Sai2.close();
}
