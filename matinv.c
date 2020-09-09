// cf. https://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
#include <stdio.h>
#include <lapacke.h>

int main() {
  int N = 3;
  int NN = 9;
  //double M[3*3] = { // ColMajor
  //  1, 2,  3,
  //  4, 5,  6,
  //  7, 8, 19};
  //double M[3*3] = { // RowMajor
  //  1,2,3, 4,5,6, 7,8,19};
  double M[3*3] = { // RowMajor
    1,4,7, 2,5,8, 3,6,19};
  int pivotArray[3]; //since our matrix has three rows
  int errorHandler;
  double lapackWorkspace[9];

  // dgetrf(M,N,A,LDA,IPIV,INFO) means invert LDA columns of an M by N matrix
  // called A, sending the pivot indices to IPIV, and spitting error information
  // to INFO. also don't forget (like I did) that when you pass a two-dimensional
  // array to a function you need to specify the number of "rows"
  dgetrf_(&N, &N, M, &N, pivotArray, &errorHandler);
  printf ("dgetrf eh, %d, should be zero\n", errorHandler);

  dgetri_(&N, M, &N, pivotArray, lapackWorkspace, &NN, &errorHandler);
  printf ("dgetri eh, %d, should be zero\n", errorHandler);

  for (size_t row = 0; row < N; ++row) {
    for (size_t col = 0; col < N; ++col) {
      printf ("%+.4f", M[row+N*col]);
      if (N-1 != col) printf (", ");
    }
    if (N-1 != row) printf ("\n");
  }
  printf ("\n");
  return 0;
}
