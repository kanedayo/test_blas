// cf. https://www.r-ccs.riken.jp/wp-content/uploads/2019/05/nakata190523.pdf
#include <stdio.h>
#include <cblas.h>
//Matlab/Octave format
void printmat(int N, int M, double *A, int LDA) {
  double mtmp;
  printf("[ ");
  for (int i = 0; i < N; i++) {
    printf("[ ");
    for (int j = 0; j < M; j++) {
      mtmp = A[i + j * LDA];
      printf("%5.2e", mtmp);
      if (j < M - 1) printf(", ");
    } if (i < N - 1) printf("]; ");
    else printf("] ");
  } printf("]");
}
int main()
{
  int n = 3; double alpha, beta;
  double *A = new double[n*n];
  double *B = new double[n*n];
  double *C = new double[n*n];
  A[0+0*n]=1; A[0+1*n]= 8; A[0+2*n]= 3;
  A[1+0*n]=2; A[1+1*n]=10; A[1+2*n]= 8;
  A[2+0*n]=9; A[2+1*n]=-5; A[2+2*n]=-1;
  B[0+0*n]= 9; B[0+1*n]= 8; B[0+2*n]=3;
  B[1+0*n]= 3; B[1+1*n]=11; B[1+2*n]=2.3;
  B[2+0*n]=-8; B[2+1*n]= 6; B[2+2*n]=1;
  C[0+0*n]=3; C[0+1*n]=3; C[0+2*n]=1.2;
  C[1+0*n]=8; C[1+1*n]=4; C[1+2*n]=8;
  C[2+0*n]=6; C[2+1*n]=1; C[2+2*n]=-2;
  printf("# dgemm demo...\n");
  printf("A =");printmat(n,n,A,n);printf("\n");
  printf("B =");printmat(n,n,B,n);printf("\n");
  printf("C =");printmat(n,n,C,n);printf("\n");
  alpha = 3.0; beta = -2.0;

  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
      n, n, n, alpha, A, n, B, n, beta, C, n);
  printf("alpha = %5.3e\n", alpha);
  printf("beta = %5.3e\n", beta);
  printf("ans="); printmat(n,n,C,n);
  printf("\n");
  printf("#check by Matlab/Octave by:\n");
  printf("alpha * A * B + beta * C =\n");
  delete[]C; delete[]B; delete[]A;
}
