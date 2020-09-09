#include <stdio.h>
#include <cblas.h>

void DebugPrint(OPENBLAS_CONST enum CBLAS_ORDER Order,
  openblas_complex_float*X, blasint M,blasint N){ // ldX==N
  for(int m=0; m<M; m++){
    for(int n=0; n<N; n++){
        openblas_complex_float x;
      if(Order==CblasRowMajor) x = X[m*((int)N)+n];
      else                     x = X[m+((int)M)*n];
      printf("[%6.1f +", openblas_complex_float_real(x)); // CblasRowMajor
      printf(" %6.1fj]", openblas_complex_float_imag(x)); // CblasRowMajor
    };printf("\n");
  };printf("\n");
}
int main(void){
  /*
  //xhat = pinv(H'*H+sig*I)*H'*y
  P=4; L=2; sig=0;
  x=(1:L).';
  H = reshape(1:P*L,P,L);H=H+1j*H;

  y = H*x;
  I = eye(L);
  hhh = H'*H + sig*I;
  inv = pinv(hhh);
  csi0 = 1./real(diag(inv));
  hhy  = H'*y;
  xhat0= inv * hhy;

  HHH = H'*H + sig*I;           % CHERK % C = A' * A + beta * C
  GH=H/HHH;                     % CHESV % X*A = B % (GH)*(HHH) = (H)
  %csi = 1./diag(GH'*GH);       % CGEMV % C = A' * B
  csi = 1./sum(conj(GH).*GH).'; % CDOTC % c = a' * b
  xhat = GH' * y ;              % CGEMV % C = A' * B

  [csi0  csi ]
  [xhat0 xhat]

  */

  //typedef struct { float real, imag; } openblas_complex_float;
  blasint ldH, incX,incY;
  blasint P=4, L=2;
  openblas_complex_float Const1={1,0}, Const0={0,0};
  openblas_complex_float sigma={0,0}; // nVar
  openblas_complex_float alpha, beta;
  // H[PxL] , X[Lx1]
  openblas_complex_float H[4*2] = { // CblasRowMajor
    {1,1}, {5,5},
    {2,2}, {6,6},
    {3,3}, {7,7},
    {4,4}, {8,8}}; ldH=2;
  openblas_complex_float X[2*1] = { // CblasRowMajor
    {1,0},
    {2,0}}; incX = 1;
  openblas_complex_float Y[4*1] = {0}; incY = 1; // CblasRowMajor

  // Y = H * X
  // Y = alpha * A * X + beta * Y // A[m,n]
  // Y = (1.0) * H * X + (0.0)* Y
  cblas_cgemv(CblasRowMajor, CblasNoTrans, P, L,
      &Const1, H, ldH,  X, incX,  &Const0, Y, incY );
  //DebugPrint(CblasRowMajor, Y, P, 1);


  // HHH = H'*H + sig*I;
  // C   = alpha * A * B + beta * C // C[m,n]
  // HHH = (1.0) * H'* H + sigma* I
  openblas_complex_float E[2*2]={1,0, 0,1}; blasint ldE=2;
  cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, L,L,P,
      &Const1, H,ldH, H,ldH, &sigma, E,ldE );
  //DebugPrint(CblasRowMajor, E, L, L);

  // inv = pinv(HHH);
  // C   = alpha * A * B + beta * C // C[m,n]
  // HHH = (1.0) * H'* H + sigma* I
  cblas_ctrsm(CblasRowMajor, CblasLeft,
      OPENBLAS_CONST enum CBLAS_SIDE Side,
      OPENBLAS_CONST enum CBLAS_UPLO Uplo,
      OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,

      OPENBLAS_CONST enum CBLAS_DIAG Diag,
      OPENBLAS_CONST blasint M,
      OPENBLAS_CONST blasint N,
      OPENBLAS_CONST void *alpha,
      OPENBLAS_CONST void *A,
      OPENBLAS_CONST blasint lda,
      void *B,
      OPENBLAS_CONST blasint ldb);
//typedef enum CBLAS_UPLO      {CblasUpper=121, CblasLower=122} CBLAS_UPLO;
//typedef enum CBLAS_DIAG      {CblasNonUnit=131, CblasUnit=132} CBLAS_DIAG;
//typedef enum CBLAS_SIDE      {CblasLeft=141, CblasRight=142} CBLAS_SIDE;

  //
  //
  // csi = 1/real(diag(inv));
  // hhy = H'*y;
  // xhat= inv * hhy

  return 0;
}
