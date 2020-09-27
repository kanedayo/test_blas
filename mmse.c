#include <stdio.h>
#include <blas.h>
#include <lapack.h>
#include <stdlib.h>

#define P_MAX 8
#define L_MAX 8
#define UPLO  "U"

void DebugPrint(complex *X, int M, int N){ // ldX==N
  for(int m=0; m<M; m++){
    for(int n=0; n<N; n++){
        complex x = X[m+(M)*n];
#if 0
      printf("[%6.1f +", x.r); // real
      printf(" %6.1fi]", x.i); // imag
#else
      printf("[ %-+11.4e,", x.r); // real
      printf(" %-+11.4ei ]", x.i); // imag
#endif
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

  HHH = H'*H + sig*I            % CHERK % C = A' * A + beta * C
  %GH=H/HHH                      % CHESV(?)% X*A = B % (GH)*(HHH) = (H)
  HU=chol(HHH)                  % CPOTRF
  GH=(H/HU)/HU'                 % CTRSM % (GH*HU') * HU = H
  %csi = 1./diag(GH'*GH)         % CGEMV % C = A' * B
  csi = 1./sum(conj(GH).*GH).'  % CDOTC % c = a' * b
  xhat = GH' * Y                % CGEMV % C = A' * B

  */

  const ptrdiff_t const0=0, const1=1, inc1 = 1; // Const
  const complex Const0[1]={{0,0}}, Const1[1]={{1,0}}; // Complex-Const

  ptrdiff_t P=4, L=2;
  ptrdiff_t ldH;
  complex sigma={0,0}; // nVar
  complex H[P_MAX*L_MAX]={0}; ldH=P;
  complex X[L_MAX] = {0};
  {
    // H[PxL] , X[Lx1]
#if 1
    //complex Hrm[4*2] = { // RowMajor
    //  {1,1}, {5,5},
    //  {2,2}, {6,6},
    //  {3,3}, {7,7},
    //  {4,4}, {8,8}};// ldH=2;
    complex Hcm[4*2] = { // ColMajor
      {1,1}, {2,2}, {3,3}, {4,4},
      {5,5}, {6,6}, {7,7}, {8,8}};// ldH=4;
    ptrdiff_t PL = P*L;
    ccopy( &PL, (float*)Hcm,&inc1, (float*)H,&inc1 ); ldH=4;

    complex Xv[2*1] = { // Vector
      {1,0},
      {2,0}};
    ccopy( &L, (float*)Xv,&inc1, (float*)X,&inc1);
#else

#define C(re,im) ((complex){re,im})
    //H[0]=C(1,1); H[4]=C(5,5);
    //H[1]=C(2,2); H[5]=C(6,6);
    //H[2]=C(3,3); H[6]=C(7,7);
    //H[3]=C(4,4); H[7]=C(8,8);
    //X[0]=C(1,1);
    //X[1]=C(2,2);
    //X[2]=C(3,3);
    //X[3]=C(4,4);

    for(int p=0; p<P_MAX; p++){
      for(int l=0; l<L_MAX; l++){
        H[P_MAX*l+p]=C(rand(),rand());
        //H[P_MAX*l+p]=C(100*rand(),100*rand());
      }
    }
    for(int l=0; l<L_MAX; l++){
      X[l]=C((l+1),0);
      //X[l]=C((l*100+100 +l*10+10 +l+1),0);
    }
#endif
  }

  complex Y[P_MAX*1] = {0};
  { // Y = H * X
    // Y = (1.0) * H * X + (0.0)* Y
    // C = alpha*A*B + beta*C // _GEMV : C[m,n]
    cgemv("N", &P, &L,     // Trans, M, N
        (float*)Const1,    // alpha
        (float*)H, &ldH,   // op(A)[M,N]
        (float*)X, &inc1,  // B
        (float*)Const0,    // beta
        (float*)&Y, &inc1);// C
    printf("DebugPrint:Y\n");
    DebugPrint( Y, P, 1 );
  }

  complex HHH[L_MAX*L_MAX]={0}; ptrdiff_t ldHHH=L;
  { // HHH = H'*H + sig*I;
    // HHH = (1.0) * H'* H + sigma* I
    // C   = alpha * A'* B + beta * C // _GEMM : C[m,n]
    // C   = alpha * A'* A + beta * C // _HERK : C[m,n]
    complex E[L_MAX*L_MAX]={{1,0},{0}, {0},{1,0}}; ptrdiff_t ldE=L;
#if 0
    cgemm("C", "N", &L,&L,&P, // Trans, Trans, M, N, K
        (float*)Const1,  // alpha
        (float*)H, &ldH, // op(A)[M,K]
        (float*)H, &ldH, // op(B)[K,N]
        (float*)&sigma,  // beta
        (float*)E, &ldE);// C[M,N]
#else
    cherk(UPLO, "C", &L,&P, // UpLo, Trans, N, K
        (float*)Const1,  // alpha
        (float*)H, &ldH, // op(A)
        (float*)&sigma,  // beta
        (float*)E, &ldE);// C
#endif
    ptrdiff_t len = L*L;
    ccopy( &len, (float*)E,&inc1, (float*)HHH,&inc1 );
    printf("DebugPrint:HHH\n");
    DebugPrint( HHH, L, L );
  }

  complex HU[L_MAX*L_MAX]={0}; ptrdiff_t ldHU=L;
  { // HU = chol(HHH);
    // HHH=HU'*HU
    // A = L*L' = U'*U // _POTRF : A[n,n]
    ptrdiff_t len = ldHHH*ldHHH;
    ccopy( &len, (float*)HHH,&inc1, (float*)HU,&inc1 );
    ptrdiff_t info;

    cpotrf(UPLO, &L,      // UpLo, N
        (float*)HU,&ldHHH,// A[N,N] -> AL or AU
        &info);
    printf("DebugPrint:HU:cpotrf.info=%d\n",(int)info);
    DebugPrint( HU, L, L );
  }

  complex GH[P_MAX*L_MAX]={0}; ptrdiff_t ldGH=P;
  { // GH=H/HU/HU'
    // (GH*HU')*HU = H % TRSM
    ptrdiff_t len = P*L;
    ccopy( &len, (float*)H,&inc1, (float*)GH,&inc1 );

    char *TRANS0 =(UPLO[0]=='U')? "N":"C";
    char *TRANS1 =(UPLO[0]=='U')? "C":"N";

    ctrsm("R",UPLO,TRANS0,"N", &P,&L,
        (float*)Const1,
        (float*)HU,&L,
        (float*)GH,&P);
    printf("DebugPrint:H/HU\n");
    DebugPrint( GH, P, L );

    ctrsm("R",UPLO,TRANS1,"N", &P,&L,
        (float*)Const1,
        (float*)HU,&L,
        (float*)GH,&P);
    printf("DebugPrint:GH=H/HU/HU'\n");
    DebugPrint( GH, P, L );
  }

  complex Xhat[L_MAX*1]={0};
  { // Xhat=GH'*Y
    // Xhat=(1.0)*A'*B + (0.0)*Xhat
    // C   = alpha*A*B + beta*C // _GEMV : C[m,n]

    cgemv("C", &P, &L,    // Trans, M, N
        (float*)Const1,   // alpha
        (float*)GH, &ldGH,// A
        (float*)Y, &inc1, // B
        (float*)Const0,   // beta
        (float*)Xhat, &inc1 );
    printf("DebugPrint:Xhat\n");
    DebugPrint( Xhat, L, 1 );
  }

  complex CSI[L_MAX*1]={0};
  { // CSI = 1./sum(conj(GH).*GH).';
    // DOT[i] = GH(:,i)' * GH(:,i) % CDOTC % c = a' * b
    // CSI[i] = 1/DOT

    for(int i=0; i<L; i++){
#if 0
      complex dot = cdotc( &P, (float*)&(GH[P*i]), &inc1, (float*)&(GH[P*i]), &inc1 );
#else
      complex dot = {0,0};
      complex*ptr = GH+(P*i);
      for(int p=0; p<P; p++,ptr++){
        float re = ptr->r;
        float im = ptr->i;
        dot.r += (re*re)+(im*im);
      }
#endif
      CSI[i].r = 1.0/dot.r;
      CSI[i].i = 0.0;
    }
    printf("DebugPrint:CSI\n");
    DebugPrint( CSI, L, 1 );

    printf("CSI=\n 3.6782\n 21.3333\n");

  }
  return 0;
}
/*
H=reshape(1:8,4,2); H=H+1j*H;
X=[1;2];
Y=H*X;

HHH=H'*H;
[HL,HU]=lu  (HHH); % HHH=HL *HU
[   HU]=chol(HHH); % HHH=HU'*HU
// _POTRF+_POTRS(LL), _SYTRF+_SYTRS(LDL),
GH=(H/HU/HL );
GH=(H/HU/HU');
CSI = 1./diag(     GH' *GH)  ; % CGEMV % C = A' * B
CSI = 1./ sum(conj(GH).*GH).'; % CDOTC % c = a' * b
xhat=GH'*Y;


% xhat
% = H\Y
% = HHH\H'*Y
% = (H/HHH)'*Y
% = (H/(HL*HU))'*Y
% = (H/HU/HL)'*Y

*/
