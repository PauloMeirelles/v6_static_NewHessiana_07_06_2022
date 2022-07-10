#ALGORITMO MEF POSICIONAL ELASTICO (ÚNICO PROCESSADOR NA CONSTRUÇÃO E SOLVER)

~~~html
# Funcinamento com as seguintes diretivas:
> Alteração do caminho do arquivo de leitura .txt

> Deletar arquivos da pasta build com "rm -r *" e "s" para remover bibliotecas baixadas na construção no "ccmake .."

> Ordenação dos nós segue:

   /|    2    
  /_|   76   
 /__|  895  
/___| 0341
~~~

#OBSEVAÇÃO: IMPLEMENTADO vtu, incremento de passo de força, nova rotina de construção matriz Hessiana tipo Petsc


**Solver em ajustes MUMPS(LU e Cholesky), novo gera .vtu e CORREÇÃO DO PROGRAMA GLOBAL "**





