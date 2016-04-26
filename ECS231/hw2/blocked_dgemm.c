#include <iostream>
#include <cstdlib>

// Get wall clock time in seconds
inline double get_time() {
  struct timespec tv;
  clock_gettime(CLOCK_MONOTONIC, &tv);
  return tv.tv_sec + tv.tv_nsec * 1e-9;
}

using namespace std;


/*
  In case you're wondering, dgemm stands for Double-precision, GEneral
  Matrix-Matrix multiplication.
*/

/*
  A is M-by-K
  B is K-by-N
  C is M-by-N

  lda is the leading dimension of the matrix (the M of square_dgemm).
*/

void
basic_dgemm (const unsigned lda,
             const unsigned M, const unsigned N, const unsigned K,
             const double *A, const double *B, double *C)
{
  unsigned i, j, k;

  /*
    To optimize this, think about loop unrolling and software
    pipelining.  Hint:  For the majority of the matmuls, you
    know exactly how many iterations there are (the block size)...
  */

  for (i = 0; i < M; ++i) {
       const double *Ai_ = A + i;
       for (j = 0; j < N; ++j) {
            const double *B_j = B + j*lda;

            double cij = *(C + j*lda + i);

            for (k = 0; k < K; ++k) {
                 cij += *(Ai_ + k*lda) * *(B_j + k);
            }

            *(C + j*lda + i) = cij;
       }
  }
}

/* You'll definitely change this... */
#define BLOCK_SIZE ((unsigned) 32)

void
do_block (const unsigned lda,
          const double *A, const double *B, double *C,
          const unsigned i, const unsigned j, const unsigned k)
{
     /*
       Remember that you need to deal with the fringes in each
       dimension.

       If the matrix is 7x7 and the blocks are 3x3, you'll have 1x3,
       3x1, and 1x1 fringe blocks.

             xxxoooX
             xxxoooX
             xxxoooX
             oooxxxO
             oooxxxO
             oooxxxO
             XXXOOOX

       You won't get this to go fast until you figure out a `better'
       way to handle the fringe blocks.  The better way will be more
       machine-efficient, but very programmer-inefficient.
     */
     const unsigned M = (i+BLOCK_SIZE > lda? lda-i : BLOCK_SIZE);
     const unsigned N = (j+BLOCK_SIZE > lda? lda-j : BLOCK_SIZE);
     const unsigned K = (k+BLOCK_SIZE > lda? lda-k : BLOCK_SIZE);

     basic_dgemm (lda, M, N, K,
                  A + i + k*lda, B + k + j*lda, C + i + j*lda);
}

void
square_dgemm (const unsigned M, 
              const double *A, const double *B, double *C)
{
     const unsigned n_blocks = M / BLOCK_SIZE + (M%BLOCK_SIZE? 1 : 0);
     unsigned bi, bj, bk;

     for (bi = 0; bi < n_blocks; ++bi) {
          const unsigned i = bi * BLOCK_SIZE;
          
          for (bj = 0; bj < n_blocks; ++bj) {
               const unsigned j = bj * BLOCK_SIZE;

               for (bk = 0; bk < n_blocks; ++bk) {
                    const unsigned k = bk * BLOCK_SIZE;
                    
                    do_block (M, A, B, C, i, j, k);
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
