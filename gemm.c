// cf. https://qiita.com/nek0log/items/a90c337ea50a11dd71d4
#include <stdio.h>
#include <cblas.h>

void print_matrix(double *data,int row,int col){
    for (int i = 0; i < row; ++i){
        for (int j = 0; j < col; ++j){
            printf("%6.1f", data[i*col + j]);
        }
        printf("\n");
    }
}

int main(void){

    double alpha, beta;
    int M = 3, N = 4, K = 2;
    // A[MxK] , B[KxN] , C[MxN]
    double A[3*2] = {
      1, 2,
      3, 4,
      5, 6};
    double B[2*4] = {
      1, 2, 3, 4,
      5, 6, 7, 8};
    double C[3*4];

    // C = alpha * A * B + beta * C
    alpha=1; beta=0; // C = A*B
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K,
    alpha,  A,K,  B,N,  beta,  C,N);

    print_matrix(C,M,N);

    // C = alpha * A * B + beta * C
    alpha=1; beta=1; // C = A*B + C
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K,
    alpha,  A,K,  B,N,  beta,  C,N);

    print_matrix(C,M,N);

    return 0;
}

/*========================
$ ./gemm.exe
  11.0  14.0  17.0  20.0
  23.0  30.0  37.0  44.0
  35.0  46.0  57.0  68.0
  22.0  28.0  34.0  40.0
  46.0  60.0  74.0  88.0
  70.0  92.0 114.0 136.0
 */
/*========================
$ octave -q --no-gui
octave:1> A=[1 2; 3 4; 5 6]; B=[1 2 3 4; 5 6 7 8];
octave:2> C=A*B, C=A*B+C
C =
   11   14   17   20
   23   30   37   44
   35   46   57   68

C =
    22    28    34    40
    46    60    74    88
    70    92   114   136
 */
