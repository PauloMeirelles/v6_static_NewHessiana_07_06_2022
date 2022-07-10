#include "MumpsSolve.h"


// ===== SOLVER EM PROJETO (NAO UTILIZAR) =========

// ===== SOLVER EM PROJETO (NAO UTILIZAR) =========

// ===== SOLVER EM PROJETO (NAO UTILIZAR) =========

// ===== SOLVER EM PROJETO (NAO UTILIZAR) =========

// ===== SOLVER EM PROJETO (NAO UTILIZAR) =========

// ===== SOLVER EM PROJETO (NAO UTILIZAR) =========

// ===== SOLVER EM PROJETO (NAO UTILIZAR) =========



#if defined(PETSC_HAVE_MUMPS)
 /* Subroutine contributed by Varun Hiremath */
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

V1D MUMPSsolver(int argc, char *argv[], V2D& Ma, V1D& Vb)
{
    // ----------------- CRIAÇÃO DAS VARIAVEIS -------------------
    // ===========================================================
    V1D  Vx;
    KSP  ksp;
    PC    pc; 
    Vec x, b, VetorB;
    Mat    A,F;
    PetscInt  N1,M1;
    //PetscInt  m1,n1;
    PetscReal      norm;     /* norm of solution error */
    PetscInt       i,j,Ii,J,Istart,Iend,m = 8,n = 7,its;
    PetscBool      flg=PETSC_FALSE,flg_ilu=PETSC_FALSE,flg_ch=PETSC_FALSE;
    //PetscInt       i,j,m = 2,n,Ii,J; //verificar
    N1= Ma[0].size();
    M1=N1/8;
    Vx.resize(N1,0.0);

 #if defined(PETSC_HAVE_MUMPS)
    PetscBool      flg_mumps=PETSC_FALSE,flg_mumps_ch=PETSC_FALSE;
 #endif
 #if defined(PETSC_HAVE_SUPERLU) || defined(PETSC_HAVE_SUPERLU_DIST)
    PetscBool      flg_superlu=PETSC_FALSE;
 #endif
 #if defined(PETSC_HAVE_STRUMPACK)
    PetscBool      flg_strumpack=PETSC_FALSE;
 #endif
    PetscScalar    v;
    PetscMPIInt    rank,size;
  #if defined(PETSC_USE_LOG)
   PetscLogStage  stage;
 #endif
    
    PetscInitialize(&argc, &argv, NULL, "RESOLVENDO UM SISTEMA ESPARÇO\n");
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);
    PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);

    // ----------------- CRIAÇÃO DA MATRIZ     -------------------
    // =========================================================== 
    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N1,N1);
    MatSetFromOptions(A);
    MatMPIAIJSetPreallocation(A,M1,NULL,M1,NULL);
    MatSeqAIJSetPreallocation(A,M1,NULL);
    MatSetUp(A); 

    MatGetOwnershipRange(A,&Istart,&Iend);

    PetscLogStageRegister("Assembly", &stage);
    PetscLogStagePush(stage);

 
        for (int i=0; i<N1; i++){
            for (int j=0; j<N1; j++){    
                const PetscInt* p_i=&i;
                const PetscInt* p_j=&j;
                const PetscScalar* p_Ma=&Ma[i][j];     
                MatSetValues(A, 1, p_i, 1, p_j, p_Ma, INSERT_VALUES); //TESTAR AINDA
            }
        }
    

    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    PetscLogStagePop();

    MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);


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


    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,A,A);
/*
     Example of how to use external package MUMPS
     Note: runtime options
           '-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_7 3 -mat_mumps_icntl_1 0.0'
           are equivalent to these procedural calls
   */
 #if defined(PETSC_HAVE_MUMPS)
   flg_mumps    = PETSC_FALSE;
   flg_mumps_ch = PETSC_FALSE;
    PetscOptionsGetBool(NULL,NULL,"-use_mumps_lu",&flg_mumps,NULL);
    PetscOptionsGetBool(NULL,NULL,"-use_mumps_ch",&flg_mumps_ch,NULL);
    if (flg_mumps || flg_mumps_ch) {
     KSPSetType(ksp,KSPPREONLY);
     PetscInt  ival,icntl;
     PetscReal val;
     KSPGetPC(ksp,&pc);
     if (flg_mumps) {
       PCSetType(pc,PCLU);
     } else if (flg_mumps_ch) {
       MatSetOption(A,MAT_SPD,PETSC_TRUE); /* set MUMPS id%SYM=1 */
       PCSetType(pc,PCCHOLESKY);
     }
     PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
     PCFactorSetUpMatSolverType(pc); /* call MatGetFactor() to create F */
     PCFactorGetMatrix(pc,&F);

     if (flg_mumps) {
       /* Get memory estimates from MUMPS' MatLUFactorSymbolic(), e.g. INFOG(16), INFOG(17).
          KSPSetUp() below will do nothing inside MatLUFactorSymbolic() */
       MatFactorInfo info;
       MatLUFactorSymbolic(F,A,NULL,NULL,&info);
       flg = PETSC_FALSE;
       PetscOptionsGetBool(NULL,NULL,"-print_mumps_memory",&flg,NULL);
       if (flg) {
         printMumpsMemoryInfo(F);
       }
     }

     /* sequential ordering */
     icntl = 7; ival = 2;
     MatMumpsSetIcntl(F,icntl,ival);

     /* threshold for row pivot detection */
     MatMumpsSetIcntl(F,24,1);
     icntl = 3; val = 1.e-6;
     MatMumpsSetCntl(F,icntl,val);

     /* compute determinant of A */
     MatMumpsSetIcntl(F,33,1);
   }
 #endif

     /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                       Solve the linear system
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   KSPSolve(ksp,b,x);

   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                       Check solution and clean up
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

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