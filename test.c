#include <stdio.h>
#include <cblas.h>

void DebugPrint(OPENBLAS_CONST enum CBLAS_ORDER Order,
  double*X, blasint M,blasint N){ // ldX==N
  for(int m=0; m<M; m++){
    for(int n=0; n<N; n++){
      if(Order==CblasRowMajor){
        printf("%6.1lf", X[m*((int)N)+n]); // CblasRowMajor
      } else {
        printf("%6.1lf", X[m+((int)M)*n]); // CblasColMajor
      }
    };printf("\n");
  };printf("\n");
}
int main(void){
  blasint ldA,ldB,ldC;
  blasint M=3, N=4, K=2;
  double alpha, beta;
    // A[MxK] , B[KxN] , C[MxN]
    double A[3*2] = { // CblasRowMajor
      1, 2,
      1, 2,
      1, 2}; ldA=2;
    double B[2*4] = { // CblasRowMajor
      1, 2, 3, 4,
      1, 2, 3, 4}; ldB = 4;
    double C[3*4] = {0}; ldC = 4; // CblasRowMajor

  // C = alpha * A * B + beta * C

  // CblasRowMajor
  alpha=1; beta=1; // C = A*B + C
  ldC = N; // CblasRowMajor
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M,N,K,
      alpha, A,ldA, B,ldB, beta, C,ldC );

  DebugPrint(CblasRowMajor,C,M,N);

  // Transpose : C(CblasRowMajor) -> D(CblasColMajor)
  alpha=1; beta=0; // C = Tr(A)*E
  double E[4*4]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1}; blasint ldE=4;
  double D[3*4];
  ldC = M; // CblasColMajor
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, M,N,N,
      alpha, C,N, E,ldE, beta, D,ldC );

  DebugPrint(CblasColMajor,D,M,N);

  // CblasColMajor
  ldC = M; // CblasColMajor
  alpha=1; beta=1; // C = A*B + C
  cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, M,N,K,
      alpha, A,ldA, B,ldB, beta, D,ldC );

  DebugPrint(CblasColMajor,D,M,N);

  return 0;
}

/*
$ make test && ./test
   3.0   6.0   9.0  12.0
   3.0   6.0   9.0  12.0
   3.0   6.0   9.0  12.0

   3.0   6.0   9.0  12.0
   3.0   6.0   9.0  12.0
   3.0   6.0   9.0  12.0

   6.0  12.0  18.0  24.0
   6.0  12.0  18.0  24.0
   6.0  12.0  18.0  24.0
 */
