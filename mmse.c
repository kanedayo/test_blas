#include <stdio.h>
#include <blas.h>
#include <lapack.h>

void DebugPrint(complex *X, int M, int N){ // ldX==N
  for(int m=0; m<M; m++){
    for(int n=0; n<N; n++){
        complex x = X[m+(M)*n];
      printf("[%6.1f +", x.r); // real
      printf(" %6.1fj]", x.i); // imag
    };printf("\n");
  };printf("\n");
}
int main(void){
  /*
  //xhat = pinv(H'*H+sig*I)*H'*y
  P=4; L=2; sig=0;
  X=(1:L).';
  H = reshape(1:P*L,P,L);H=H+1j*H;

  Y = H*X; I = eye(L);

  // HHH = H'*H + sig*I;
  // inv = pinv(HHH);
  // csi0 = 1./real(diag(inv));
  // hhy  = H'*y;
  // xhat0= inv * hhy;

  HHH = H'*H + sig*I;           % CHERK % C = A' * A + beta * C
  GH=H/HHH;                     % CHESV % X*A = B % (GH)*(HHH) = (H)
  %csi = 1./diag(GH'*GH);       % CGEMV % C = A' * B
  csi = 1./sum(conj(GH).*GH).'; % CDOTC % c = a' * b
  xhat = GH' * y ;              % CGEMV % C = A' * B

  [csi0  csi ]
  [xhat0 xhat]

  */

  ptrdiff_t ldH, incX,incY;
  ptrdiff_t P=4, L=2;
  complex Const1={1,0}, Const0={0,0};
  complex sigma={0,0}; // nVar
  //complex alpha, beta;
  complex H[4*2]={0}; ldH=0;
  complex X[2] = {0}; incX = 1;
  complex Y[4] = {0}; incY = 1;
  {
    // H[PxL] , X[Lx1]
    //complex Hrm[4*2] = { // RowMajor
    //  {1,1}, {5,5},
    //  {2,2}, {6,6},
    //  {3,3}, {7,7},
    //  {4,4}, {8,8}};// ldH=2;
    complex Hcm[4*2] = { // ColMajor
      {1,1}, {2,2}, {3,3}, {4,4},
      {5,5}, {6,6}, {7,7}, {8,8}};// ldH=4;
    ptrdiff_t PL = P*L;
    ptrdiff_t incH = 1;
    ccopy( &PL, (float*)Hcm,&incH, (float*)H,&incH ); ldH=4; // ColMajor

    complex Xv[2*1] = { // Vector
      {1,0},
      {2,0}};// incX = 1;
    //X = Xv; incX = 1;
    ccopy( &L, (float*)Xv,&incX, (float*)X,&incX);

    complex Yv[4*1] = {0};// incY = 1; // Vector
    //Y = Yv; incY = 1;
    ccopy( &P, (float*)Yv,&incY, (float*)Y,&incY);
  }

  { // Y = H * X
    // Y = alpha * A * X + beta * Y // A[m,n]
    // Y = (1.0) * H * X + (0.0)* Y
    cgemv("N", &P, &L,
        (float*)&Const1,
        (float*)H, &ldH,
        (float*)X, &incX,
        (float*)&Const0,
        (float*)&Y, &incY );
    DebugPrint( Y, P, 1 );
  }

  complex HHH[2*2]={0}; ptrdiff_t ldHHH=2;
  { // HHH = H'*H + sig*I;
    // C   = alpha * A * B + beta * C // C[m,n]
    // HHH = (1.0) * H'* H + sigma* I
    complex E[2*2]={{1,0},{0}, {0},{1,0}}; ptrdiff_t ldE=2;
#if 1
    cgemm("CT", "N", &L,&L,&P,
        (float*)&Const1,
        (float*)H, &ldH,
        (float*)H, &ldH,
        (float*)&sigma,
        (float*)E, &ldE );
#else
    cherk("L", "CT", &L,&P,
        (float*)&Const1,
        (float*)H, &ldH,
        (float*)&sigma,
        (float*)E, &ldE );
#endif
    ptrdiff_t len = 2*2;
    ptrdiff_t inc = 1;
    ccopy( &len, (float*)E,&inc, (float*)HHH,&inc );
    DebugPrint( E, L, L );
  }

  complex GH[2*2]={0}; ptrdiff_t ldGH=2;
  { // GH=H/HHH;                     % CHESV % X*A = B % (GH)*(HHH) = (H)
    // CHESV, CPOSV, TRSM , _SYTRF+_SYTRS

    ptrdiff_t ipiv[2*2]={0};
    float work[16*16];
    ptrdiff_t lwork = sizeof(work);
    ptrdiff_t info;
    csytrf(
        "U", &L,
        (float*)HHH, &L, // HHH => L*L'
        ipiv,
        work, &lwork,
        &info);
    printf("DebugPrint:L\n");
    DebugPrint( HHH, L, L );


    {
      complex TMP[2*2]={0}; ptrdiff_t ldTMP=2;
      ptrdiff_t len = 2*2;
      ptrdiff_t inc = 1;
      ccopy( &len, (float*)HHH,&inc, (float*)TMP,&inc );

      ctrmm( "R","L","H","N", &L,&L,
          (float*)&Const1, // alpha
          (float*)HHH,&ldTMP,
          (float*)TMP,&ldTMP);
      printf("DebugPrint:L*L'\n");
      DebugPrint( TMP, L, L );

    }


    csytrs(
        "U",&L, &L,
        (float*)HHH, &L,
        ipiv,
        (float*)GH, &L,
        &info);
    DebugPrint( GH, L, L );


    /*
       extern void csytrf(
       const char   *uplo,
       const ptrdiff_t *n,
       float  *a,
       const ptrdiff_t *lda,
       ptrdiff_t *ipiv,
       float  *work,
       const ptrdiff_t *lwork,
       ptrdiff_t *info
       );
       */
      /*
         extern void csytrs(
         const char   *uplo,
         const ptrdiff_t *n,
         const ptrdiff_t *nrhs,
         const float  *a,
         const ptrdiff_t *lda,
         const ptrdiff_t *ipiv,
         float  *b,
         const ptrdiff_t *ldb,
         ptrdiff_t *info
         );
         */

    /*
       extern void ctrsm(
       const char   *side,
       const char   *uplo,
       const char   *transa,
       const char   *diag,
       const ptrdiff_t *m,
       const ptrdiff_t *n,
       const float  *alpha,
       const float  *a,
       const ptrdiff_t *lda,
       float  *b,
       const ptrdiff_t *ldb
       );
       */

    /*
       extern void cposv(
       const char   *uplo,
       const ptrdiff_t *n,
       const ptrdiff_t *nrhs,
       float  *a,
       const ptrdiff_t *lda,
       float  *b,
       const ptrdiff_t *ldb,
       ptrdiff_t *info
       );
       */
  }

  //
  //
  // csi = 1/real(diag(inv));
  // hhy = H'*y;
  // xhat= inv * hhy

  return 0;
}
