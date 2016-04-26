#include <iostream>
#include <cstdlib>

// Get wall clock time in seconds
inline double get_time() {
  struct timespec tv;
  clock_gettime(CLOCK_MONOTONIC, &tv);
  return tv.tv_sec + tv.tv_nsec * 1e-9;
}

using namespace std;

// Prototype for DGEMM Fortran interface (all parameters by reference)
extern "C" {
void dgemm_(const char *transa, const char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
}

int main(int argc, const char *argv[])
{
   
  int i, j, n;
  for (n = 32; n < 1025; n=n+32)
  {
  // build test matrix A
  double *A = new double[n * n];
  for (j = 0; j < n ; j++)
    for (i = 0; i < n; i++) 
      A[j * n + i] = double(rand()) / RAND_MAX;
      
  // build test matrix B
  double *B = new double[n * n];
  for (j = 0; j < n ; j++)
    for (i = 0; i < n; i++) 
      B[j * n + i] = double(rand()) / RAND_MAX;

  // Allocate C
  double *C = new double[n * n];
  
  // Perform the matrix multiplication
  
  double alpha = 1.0, beta = 0.0;
  double t1 = get_time(); // start time
  dgemm_("N", "N", &n, &n, &n, &alpha, A, &n, B, &n, &beta, C, &n);
  double t2 = get_time();
  
  cout << n << " " << t2 - t1 << endl;
  }

  return 0;
}
