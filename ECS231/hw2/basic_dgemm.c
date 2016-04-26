#include <iostream>
#include <cstdlib>

// Get wall clock time in seconds
inline double get_time() {
  struct timespec tv;
  clock_gettime(CLOCK_MONOTONIC, &tv);
  return tv.tv_sec + tv.tv_nsec * 1e-9;
}

using namespace std;

void
square_dgemm (const unsigned M, 
              const double *A, const double *B, double *C)
{
  unsigned i, j, k;

  for (i = 0; i < M; ++i) {
    for (j = 0; j < M; ++j) {
      for (k = 0; k < M; k+=4) {
        *(C + j*M + i) += *(A + i + k*M) * *(B + j*M + k);
        *(C + j*M + i) += *(A + i + (k+1)*M) * *(B + j*M + (k+1));
        *(C + j*M + i) += *(A + i + (k+2)*M) * *(B + j*M + (k+2));
        *(C + j*M + i) += *(A + i + (k+3)*M) * *(B + j*M + (k+3));
      }
    }
  }
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
  square_dgemm(n, A, B, C);
  double t2 = get_time();

  cout <<  n  << " " << t2 - t1 << endl;
  delete[] A;
  delete[] B;
  delete[] C;
  }
  return 0;
}
