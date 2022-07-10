#include "solution.h"

/*
V1D StaticEq(int argc, char *argv[], tvEl& El, const int& ngl, const double& nFS, const double &FStep)
{
    V1D g_res; 
    std::vector<double> _g;
    std::vector<std::vector<double>> Fint;
    Fint.resize(ngl, std::vector<double>(1, 0.0));
    std::vector<std::vector<double>> Fext;
    Fext.resize(ngl, std::vector<double>(1, 0.0));
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
                g_res[El[i].r_No(a)->r_Adress(z)] = Fext[El[i].r_No(a)->r_Adress(z)][0] - Fint[El[i].r_No(a)->r_Adress(z)][0];
            }
          }
        }
    }
    return g_res;
}
*/
//SOLVER EINGEN

V1D EIGENsolver(int argc, char *argv[], V2D& Hess, V1D& _g)
{
    V1D  deltay;
    int N1;
    N1= _g.size();
    deltay.resize(N1,0.0);

    Eigen::MatrixXf MatrixEigen(N1,N1);
    for(int i=0; i< N1; i++)
    {
        for(int j=0; j<N1; j++)
        {
          MatrixEigen(i,j)=Hess[i][j];
            
        }
    }

    Eigen::VectorXf VectorEigen(N1);
    for(int i=0; i< N1; i++)
    {
        VectorEigen(i)=_g[i];
    }
    //MatrixXf m = MatrixXf::Random(n, n);
    //VectorXf b = VectorXf::Random(n);
    //auto start = sc.now();     // start timer
    Eigen::VectorXf x = MatrixEigen.lu().solve(VectorEigen);
    //auto end = sc.now();

    for(int i=0; i< N1; i++)
    {
        double p=0.0;
        p=x(i);
        deltay[i]=p;
    }
    return deltay;
}


V1D PETScSolver(int argc, char *argv[], Mat& H_static, V1D&  g_)
{
// ----------------- CRIAÇÃO DAS VARIAVEIS -------------------
    // ===========================================================
    V1D  deltay;
    KSP  ksp;
    PC    pc; 
    Vec x, b, VetorB, All;
    Mat    A;
    PetscErrorCode ierr;
    PetscInt  N1; //Istart, Iend, iteration;
    N1= g_.size();
    deltay.resize(N1,0.0);
    
    
    //VecScatter ctx;
    //PetscScalar val;
    //PetscInitialize(&argc, &argv, NULL, "RESOLVENDO UM SISTEMA ESPARÇO\n");
    
    // ----------------- CRIAÇÃO DOS VETORES   -------------------
    // ===========================================================
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, N1);
    VecSetFromOptions(x);
    VecDuplicate(x, &b);
    
    for (int i=0; i<N1; i++)
    {
        const PetscInt* p_i=&i;
        const PetscScalar* p_Vb=&g_[i];
        VecSetValues(b, 1 , p_i, p_Vb, ADD_VALUES); // ATRIBUIÇÃO DE VALORES DO VETOR Vb ENVIADO PARA O VETOR B
    }

    VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    // "b" é tipo struct (grupo de elementos de dados referentes ao vetor b)
    // std::cout << "[";
    // std::cout << *b << ", ";
    // std::cout << "]" << "\n";

    //VecView(b, PETSC_VIEWER_STDOUT_WORLD);

    VecAssemblyBegin (b); //Processo de montagem 
    VecAssemblyEnd (b);

    // ----------------- CRIAÇÃO DAS MATRIZES  -------------------
    // ===========================================================
    
    MatCreate(PETSC_COMM_WORLD, &A);
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, &A)
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N1, N1);
    MatSetFromOptions(A);
    MatSetUp(A);
    
    A=H_static;
    //MatGetOwnershipRange(A, &Istart,&Iend);
    //MatSetValues(A, 1,&dof_i, 1,&dof_j,&Ajac(2*i ,2*j),ADD_VALUES);

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    MatView(A, PETSC_VIEWER_STDOUT_WORLD);


    
    // ----------------- KSP RESOLVE -------------------
    // ===========================================================
    
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetFromOptions(ksp);
   // KSPSetTolerances(ksp, 1.e-7,1.e-10, PETSC_DEFAULT, 10000);
   // KSPGMRESSetRestart(ksp, 8);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCJACOBI);
    //KSPSetPCSide(ksp, PC_RIGHT);
    //KSPSetType(ksp, KSPLSQR);
    //KSPSetType(ksp, KSPCG);
    //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    KSPSolve(ksp, b, x);


    //VecScatterCreateToAll(x, &ctx, &All);
    //VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
    //VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);


    for(int i=0; i<N1 ;i++)
    {
        const PetscInt* p_i=&i; 
        PetscScalar vTeste;
        PetscScalar *p_vTeste=&vTeste;

        VecGetValues(x, 1, p_i, p_vTeste);
        deltay[i]=vTeste;
    }

    //VecView(x, PETSC_VIEWER_STDOUT_WORLD);

    // std::cout << "Vetor de incremento Dy" << "\n";
    // for(int h=0; h< N1; h++)
    // {
    // std::cout << deltay[h] << "\n" ;  
    // }

    //VecScatterDestroy(&ctx);
    VecDestroy(&b); //Desaloca as variaveis
    VecDestroy(&x); //Desaloca as variaveis
    MatDestroy(&A);
    KSPDestroy(&ksp);
    return deltay;
}









//================== ABRE AQUI ====================
/*
#if defined(PETSC_HAVE_MUMPS)
PetscErrorCode printMumpsMemoryInfo(Mat F)
{
     PetscInt       maxMem, sumMem;

    MatMumpsGetInfog(F,16,&maxMem);
    MatMumpsGetInfog(F,17,&sumMem);
    PetscPrintf(PETSC_COMM_WORLD, "\n MUMPS INFOG(16) :: Max memory in MB = %d", maxMem);
    PetscPrintf(PETSC_COMM_WORLD, "\n MUMPS INFOG(17) :: Sum memory in MB = %d \n", sumMem);
   return 0;
}
#endif

static PetscErrorCode MatScaleUserImpl_SeqAIJ(Mat inA,PetscScalar alpha)
{
   MatScale(inA,alpha);
   return 0;
}

extern PetscErrorCode MatScaleUserImpl(Mat,PetscScalar);

static PetscErrorCode MatScaleUserImpl_MPIAIJ(Mat A,PetscScalar aa)
{
    Mat            AA,AB;

   MatMPIAIJGetSeqAIJ(A,&AA,&AB,NULL);
   MatScaleUserImpl(AA,aa);
   MatScaleUserImpl(AB,aa);
   return 0;
}

PetscErrorCode RegisterMatScaleUserImpl(Mat A)
{
   PetscMPIInt    size;

   MPI_Comm_size(PetscObjectComm((PetscObject)A), &size);
   if (size == 1) { // SeqAIJ Matrix 
     PetscObjectComposeFunction((PetscObject)A,"MatScaleUserImpl_C",MatScaleUserImpl_SeqAIJ);
   } else { // MPIAIJ Matrix 
     Mat AA,AB;
     MatMPIAIJGetSeqAIJ(A,&AA,&AB,NULL);
     PetscObjectComposeFunction((PetscObject)A,"MatScaleUserImpl_C",MatScaleUserImpl_MPIAIJ);
     PetscObjectComposeFunction((PetscObject)AA,"MatScaleUserImpl_C",MatScaleUserImpl_SeqAIJ);
     PetscObjectComposeFunction((PetscObject)AB,"MatScaleUserImpl_C",MatScaleUserImpl_SeqAIJ);
   }
   return 0;
}


PetscErrorCode MatScaleUserImpl(Mat A,PetscScalar a)
{
    PetscErrorCode (*f)(Mat,PetscScalar);

   PetscObjectQueryFunction((PetscObject)A,"MatScaleUserImpl_C",&f);
   if (f) (*f)(A,a);
    return 0;
}

*/
/*
V1D MUMPSresolve(int argc, char *argv[], V2D& Ma, V1D& Vb)
{
    // ----------------- CRIAÇÃO DAS VARIAVEIS -------------------
    // ===========================================================
    V1D  Vx;
    KSP  ksp;
    PC    pc; 
    Vec x, b, VetorB;
    Mat    A, F;
    PetscInt  N1,n;
    //PetscInt       i,j,m = 2,n,Ii,J; //verificar
    N1= Ma[0].size();
    Vx.resize(N1,0.0);
    PetscScalar    one = 1.0;
    PetscMPIInt    rank,size;
    PetscReal      norm,tol=100.*PETSC_MACHINE_EPSILON; 
    PetscBool      flg=PETSC_FALSE,flg_ilu=PETSC_FALSE,flg_ch=PETSC_FALSE;
    PetscBool      flg_mumps=PETSC_FALSE,flg_mumps_ch=PETSC_FALSE;
    
    PetscInitialize(&argc, &argv, NULL, "RESOLVENDO UM SISTEMA ESPARÇO\n");
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);
   
    // ----------------- CRIAÇÃO DOS VETORES   -------------------
    // ===========================================================
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, N1);
    VecSetFromOptions(x);
    VecDuplicate(x, &b);
    for (int i=0; i<N1; i++)
    {
        const PetscInt* p_i=&i;
        const PetscScalar* p_Vb=&Vb[i];
        VecSetValues(b, 1 , p_i, p_Vb, ADD_VALUES) ; // ATRIBUIÇÃO DE VALORES DO VETOR Vb ENVIADO PARA O VETOR B
    }
    VecAssemblyBegin (b); //Processo de montagem 
    VecAssemblyEnd (b);

    //VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    // "b" é tipo struct (grupo de elementos de dados referentes ao vetor b)
    // std::cout << "[";
    // std::cout << *b << ", ";
    // std::cout << "]" << "\n";
  
    // ----------------- CRIAÇÃO DAS MATRIZES  -------------------
    // ===========================================================
    
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N1, N1);
    //MatSetFromOptions(A);
    MatSetType(A,MATAIJ);
    MatSetUp(A);
    //RegisterMatScaleUserImpl(A);
    
    for (int i=0; i<N1; i++)
    {
        for (int j=0; j<N1; j++)
        {    
                const PetscInt* p_i=&i;
                const PetscInt* p_j=&j;
                const PetscScalar* p_Ma=&Ma[i][j];     
                MatSetValues(A, 1, p_i, 1, p_j, p_Ma, INSERT_VALUES); //TESTAR AINDA
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    // ----------------- KSP RESOLVE -------------------
    // ===========================================================
    
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);

    #if defined(PETSC_HAVE_MUMPS)
    KSPSetType(ksp,KSPPREONLY);
    KSPGetPC(ksp,&pc);
    PCSetType(pc, PCLU);
    #endif
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
    //KSPGetPC(ksp,&pc);
    //PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
    //PCFactorSetUpMatSolverType(pc);
    //PCFactorGetMatrix(pc,&F);

    //PCSetType(pc,PCLU);
    KSPSolve(ksp, b, x);
    //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

    //PetscInt info1;
    //MatMumpsGetInfo(F,1,&info1);

    for(int i=0; i<N1 ;i++)
    {
        const PetscInt* p_i=&i; 
        PetscScalar vTeste;
        PetscScalar *p_vTeste=&vTeste;

        VecGetValues(x, 1, p_i, p_vTeste);
        Vx[i]=vTeste;
    }

    //VecView(x, PETSC_VIEWER_STDOUT_WORLD);
    
    VecDestroy(&b); //Limpa as variaveis
    VecDestroy(&x); //Limpa as variaveis
    MatDestroy(&A);
    KSPDestroy(&ksp);
    //PetscFinalize(); 
    return Vx; 
}

*/
 
//================== FECHA AQUI ====================

/*
V1D MUMPSresolve(int argc, char *argv[], V2D& Ma, V1D& Vb)
{
    // ----------------- CRIAÇÃO DAS VARIAVEIS -------------------
    // ===========================================================
    V1D  Vx;
    KSP  ksp;
    PC    pc; 
    Vec x, b, VetorB;
    Mat    A,F;
    PetscInt  N1,n;
    PetscInt  m1,n1;
    //PetscInt       i,j,m = 2,n,Ii,J; //verificar
    N1= Ma[0].size();
    Vx.resize(N1,0.0);
    PetscScalar    one = 1.0;
    PetscMPIInt    rank,size;
    PetscReal      norm,tol=100.*PETSC_MACHINE_EPSILON; 
    PetscBool      flg=PETSC_FALSE,flg_ilu=PETSC_FALSE,flg_ch=PETSC_FALSE;
    PetscBool      flg_mumps=PETSC_FALSE,flg_mumps_ch=PETSC_FALSE;
    PetscInt info1;
    //MatMumpsGetInfo(F,1,&info1); <PARA ACESSAR A FLAG INFO(!)>
    
    PetscInitialize(&argc, &argv, NULL, "RESOLVENDO UM SISTEMA ESPARÇO\n");
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);
   
    // ----------------- CRIAÇÃO DOS VETORES   -------------------
    // ===========================================================
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, N1);
    VecSetFromOptions(x);
    VecDuplicate(x, &b);
    for (int i=0; i<N1; i++)
    {
        const PetscInt* p_i=&i;
        const PetscScalar* p_Vb=&Vb[i];
        VecSetValues(b, 1 , p_i, p_Vb, ADD_VALUES) ; // ATRIBUIÇÃO DE VALORES DO VETOR Vb ENVIADO PARA O VETOR B
    }
    VecAssemblyBegin (b); //Processo de montagem 
    VecAssemblyEnd (b);

    //VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    // "b" é tipo struct (grupo de elementos de dados referentes ao vetor b)
    // std::cout << "[";
    // std::cout << *b << ", ";
    // std::cout << "]" << "\n";
  
    // ----------------- CRIAÇÃO DAS MATRIZES  -------------------
    // ===========================================================
    
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N1, N1);
    MatSetType(A,MATAIJ);
    MatSetFromOptions(A);
    MatSetUp(A);
    //RegisterMatScaleUserImpl(A);
    
    //MatSetValues(A,m1,PetscInt idxm[i],n1,PetscInt idxn[j], Ma[],ADD_VALUES);
    
    for (int i=0; i<N1; i++)
    {
        for (int j=0; j<N1; j++)
        {    
                const PetscInt* p_i=&i;
                const PetscInt* p_j=&j;
                const PetscScalar* p_Ma=&Ma[i][j];     
                MatSetValues(A, 1, p_i, 1, p_j, p_Ma, INSERT_VALUES); //TESTAR AINDA
        }
    }
    

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    // ----------------- KSP RESOLVE -------------------
    // ===========================================================
    
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    #if defined(PETSC_HAVE_MUMPS)
    KSPSetType(ksp,KSPPREONLY);
    KSPGetPC(ksp,&pc);
    PCSetType(pc, PCLU);   
    #endif
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);

    //PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
    //PCFactorSetUpMatSolverType(pc);
    PCFactorGetMatrix(pc,&A);
    KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    KSPSolve(ksp, b, x);

    MatMumpsGetInfo(A,1,&info1);


    for(int i=0; i<N1 ;i++)
    {
        const PetscInt* p_i=&i; 
        PetscScalar vTeste;
        PetscScalar *p_vTeste=&vTeste;

        VecGetValues(x, 1, p_i, p_vTeste);
        Vx[i]=vTeste;
    }

    //VecView(x, PETSC_VIEWER_STDOUT_WORLD);
    
    VecDestroy(&b); //Limpa as variaveis
    VecDestroy(&x); //Limpa as variaveis
    MatDestroy(&A);
    KSPDestroy(&ksp);
    //PetscFinalize(); 
    return Vx; 
}
*/

V2D pseudoInverse(V2D& Ma) 
{  
    double epsilon = std::numeric_limits<double>::epsilon();
    int n, m;
    n=Ma.size();
    m=Ma[1].size();
    Eigen::MatrixXd mEigen(n,m);
    for(int i=0; i< n; i++)
    {
        for(int j=0; j<m; j++)
        {
          mEigen(i,j)=Ma[i][j];
            
        }
    }
    
    Eigen::JacobiSVD< Eigen::MatrixXd > svd(mEigen ,Eigen::ComputeThinU | Eigen::ComputeThinV);  
    double tolerance = epsilon * std::max(mEigen.cols(), mEigen.rows()) *svd.singularValues().array().abs()(0); 
    Eigen::MatrixXd pinv=svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint(); 
    
    V2D invM;
    invM.resize(m,V1D (n,0.0));
    for(int i=0; i< m; i++)
    {
        for(int j=0; j<n; j++)
        {
        double p=pinv(i,j);
        invM[i][j]=p;
        }
    }

    return invM; 
}

V2D invAnyEIGEN(V2D& Ma)
{
    int n, m;
    n=Ma.size();
    m=Ma[1].size();

    Eigen::MatrixXd MatrixEigen(n,m);
    for(int i=0; i< n; i++)
    {
        for(int j=0; j<m; j++)
        {
          MatrixEigen(i,j)=Ma[i][j];
            
        }
    }

    Eigen::HouseholderQR<Eigen::MatrixXd> qr(MatrixEigen.transpose());
    Eigen::MatrixXd pinv;
    pinv.setIdentity(MatrixEigen.cols(), MatrixEigen.rows());
    pinv = qr.householderQ() * pinv;
    pinv = qr.matrixQR().topLeftCorner(MatrixEigen.rows(),MatrixEigen.rows()).triangularView<Eigen::Upper>().transpose().solve<Eigen::OnTheRight>(pinv);

    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cqr(MatrixEigen);
    pinv = cqr.pseudoInverse();

    V2D invM;
    invM.resize(m,V1D (n,0.0));
    for(int i=0; i< m; i++)
    {
      for(int j=0; j<n; j++)
      {
        double a=pinv(i,j);
        invM[i][j]=a;
      }
    }
  return invM;
}

V2D invEIGEN(V2D& Ma)
{
    
    int n, m;
    n=Ma.size();
    m=Ma[1].size();
    Eigen::MatrixXd MatrixEigen;
    MatrixEigen.resize(n,m);
    for(int i=0; i< n; i++)
    {
        for(int j=0; j<m; j++)
        {
          MatrixEigen(i,j)=Ma[i][j];
            
        }
    }
    V2D invMatrix;
    invMatrix.resize(n,V1D (m,0.0));
    Eigen::MatrixXd invMatrixEigen=MatrixEigen.inverse();
    for(int i=0; i< n; i++)
    {
        for(int j=0; j<m; j++)
        {
        double a=invMatrixEigen(i,j);
        invMatrix[i][j]=a;
        }
    }

    return invMatrix;
};

V1D KSPresolve(int argc, char *argv[], V2D& Ma, V1D& Vb)
{
    // ----------------- CRIAÇÃO DAS VARIAVEIS -------------------
    // ===========================================================
    V1D  Vx;
    KSP  ksp;
    PC    pc; 
    Vec x, b, VetorB, All;
    Mat    A;
    PetscErrorCode ierr;
    PetscInt  N1; //Istart, Iend, iteration;
    N1= Ma[0].size();
    Vx.resize(N1,0.0);
    
    //VecScatter ctx;
    //PetscScalar val;
    //PetscInitialize(&argc, &argv, NULL, "RESOLVENDO UM SISTEMA ESPARÇO\n");
    
    // ----------------- CRIAÇÃO DOS VETORES   -------------------
    // ===========================================================
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, N1);
    VecSetFromOptions(x);
    VecDuplicate(x, &b);
    
    for (int i=0; i<N1; i++)
    {
        const PetscInt* p_i=&i;
        const PetscScalar* p_Vb=&Vb[i];
        VecSetValues(b, 1 , p_i, p_Vb, ADD_VALUES); // ATRIBUIÇÃO DE VALORES DO VETOR Vb ENVIADO PARA O VETOR B
    }

    //VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    // "b" é tipo struct (grupo de elementos de dados referentes ao vetor b)
    // std::cout << "[";
    // std::cout << *b << ", ";
    // std::cout << "]" << "\n";

    VecAssemblyBegin (b); //Processo de montagem 
    VecAssemblyEnd (b);
    
    // ----------------- CRIAÇÃO DAS MATRIZES  -------------------
    // ===========================================================
    
    MatCreate(PETSC_COMM_WORLD, &A);
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, &A)
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N1, N1);
    MatSetFromOptions(A);
    MatSetUp(A);
    
    for (int i=0; i<N1; i++)
    {
        for (int j=0; j<N1; j++)
        {    
                const PetscInt* p_i=&i;
                const PetscInt* p_j=&j;
                const PetscScalar* p_Ma=&Ma[i][j];     
                MatSetValues(A, 1, p_i, 1, p_j, p_Ma, ADD_VALUES); //TESTAR AINDA
        }
    }
    
    //MatGetOwnershipRange(A, &Istart,&Iend);
    //MatSetValues(A, 1,&dof_i, 1,&dof_j,&Ajac(2*i ,2*j),ADD_VALUES);

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    //MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    

    // ----------------- KSP RESOLVE -------------------
    // ===========================================================
    
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetFromOptions(ksp);
    KSPSetTolerances(ksp, 1.e-7,1.e-10, PETSC_DEFAULT, 10000);
    KSPGMRESSetRestart(ksp, 8);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCJACOBI);
    //KSPSetPCSide(ksp, PC_RIGHT);
    //KSPSetType(ksp, KSPLSQR);
    //KSPSetType(ksp, KSPCG);
    //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    KSPSolve(ksp, b, x);


    //VecScatterCreateToAll(x, &ctx, &All);
    //VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
    //VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);


    for(int i=0; i<N1 ;i++)
    {
        const PetscInt* p_i=&i; 
        PetscScalar vTeste;
        PetscScalar *p_vTeste=&vTeste;

        VecGetValues(x, 1, p_i, p_vTeste);
        Vx[i]=vTeste;
    }

    //VecView(x, PETSC_VIEWER_STDOUT_WORLD);
    //VecScatterDestroy(&ctx);
    VecDestroy(&b); //Desaloca as variaveis
    VecDestroy(&x); //Desaloca as variaveis
    MatDestroy(&A);
    KSPDestroy(&ksp);
    return Vx;
};

