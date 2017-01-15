#include <stdio.h> 
#include <stdlib.h> 
#include <assert.h>
#include <time.h>
#define RANDOM_CONSTANT 50

/*
 * function: genMatrix( int m, int n, double *A) 
 * Usage: genMatrix( 9, 6, matrix ); 
 * Description: This function takes two integer m and n as inputs, and generate a
 * 		matrix with m by n size with random doubles, and stores it in output A
 *
 * Note: user must pass a pointer of double to this function in order for it to create
 * 	     a matrix, also caller is responsible for matrix memory allocation and memory 
 *       cleanup
 */
void genMatrix( int m, int n, double *A)
{
    // in this pointer format, we simply use address of A[0] to A[N] to represent
    // first row of data, and A[N+1] to A[2N] to represent 2nd row, ect.
    // for the purpose of this function, we just need to popular every entries
    // within this matrix with random generated double.
    double *p;
    srand48(time(0)); 
    for ( p = A; p < (A + (m*n)); p++ )
    {
        *p = drand48() * RANDOM_CONSTANT;
    }
    // by now, this function has completed its job.
}; 

/*
 * function: genVector( int n, double *b )
 * Usage:   genVectore( 6, vector ); 
 * Description: This function takes one integer, and generate a vectore with n random doubles, 
 *              and stores it in output b. 
 * Note:  User must pass a pointer of double to this function in order for it to create a matrix,
 *        also, caller is responsible for vector memory allocation and memory cleanup.
 */
void genVector( int n, double *b)
{
    double *p;
    srand48(time(0));
    for (p = b; p < ( b + n); p++)
    {
        *p = drand48() * RANDOM_CONSTANT;
    }
};
