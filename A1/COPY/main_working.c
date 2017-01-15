#include <stdio.h>
#include <string.h> 
#include "mpi.h" 
#include "genmatvec.c" 

/*
 * function: productDistributor( double *A, double *b, int n, int numP)
 * Usage: productDistributor( A, b, 10, 4);
 * Description: This function distributes the matrix A and vector b evenly based on number of
 *              available processes
 * Notes: Matrix A will be divided by column, and caller is responsible to make sure number of
 *        column is the same as number of entries of vector b so that matrix-vector mupltiplication
 *        actually makes sense.
 *Output:
 *        numofEntries: informs how many entries has been passed in a process.
 *        AA: an array that stores the column data from A that a process will be responsible for
 *        bb: an array that stores the entries of b that a process will be responsible for
 */
void productDistributor( double *A, double *b, int m, int n, int * numofEntries, double *AA, double *bb)
{
    int my_rank;            // rank of process
    int numP;               // number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numP);
    
    *numofEntries = 0; //reset number of entries.
    // assign each of entry based on size of n
    for (int i = 0; i < n; i++)
    {
        if( my_rank == i%numP ) // the assgined process which takes on the corresponding ith column and ith entry
        {
            // we need to push m entries to array AA to represent the entire column of A, or we can convert matrix A as
            // coloum based matrix
            for (int j = 0; j < m; j++)
            {
                *(AA + (*numofEntries)*m + j) = *(A + i + (j*n));
            }
            
            *(bb + *numofEntries) = *(b + i);
            *numofEntries += 1;    // to signal this process now contains one more entry
        }
    }
    
    // gmo debug print my values
    printf(" this outputs the messages for each process-- process: %i \n", my_rank);
    
    for (int i = 0; i < *numofEntries * m; i++) {
        printf(" the %i of entries for matrix is: %f \n", i, *(AA + i));
    }
    
    for (int i =0; i< *numofEntries; i++) {
        printf(" the %i of entries for vector is: %f \n", i, *(bb+i));
    }
    /* end of gmo testing */
}

/*
 * function: internalProcessMultiplication( double *AA, double *bb, int m, int numofEntries, double *result)
 * Usage: internalProcessMultiplication( AA, bb, m, numofEntries, result); 
 * Description: This function takes inputs AA: an array that contains column entries from the global matrix,
 *              bb: an array that contains entries that assigned to a process from the global vector
 *              m: number of rows from global matrx 
 *              numofEntries: number of entries that assigned to the process
 *              and computes the multiplication between sub matrix AA and sub vector bb, and stores the result 
 *              in an output array result.
 */
void internalProcessMultiplication( double *AA, double *bb, int m, int numofEntries, double *result)
{
    //Note: size of result is crucial to this computation.
    //gmo
    for (int i=0; i < numofEntries; i++)
    {
        // multiplication for each column
        for (int j = 0; j < m; j++)
        {
            *(result + i*m + j) = (*(AA + i*m + j)) * (*(bb + i));
        }
    }

    // gmo debug print my values
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    printf(" this outputs the results for each process-- process: %i \n", my_rank);
    
    for (int i = 0; i < numofEntries; i++)
    {
        for (int j=0; j < m; j++) {
            printf(" the %i of result is: %f \n", i, *(result + i*m + j));
        }
        
    }
    /* end of gmo testing */
}

int main(int argc, char* argv[])
{
    int         my_rank;       /* rank of process      */
    int         p;             /* number of processes  */
    int         source;        /* rank of sender       */
    int         dest;          /* rank of receiver     */
    int         tag = 0;       /* tag for messages     */
    char        message[100];  /* storage for message  */
    MPI_Status  status;        /* status for receive   */

    /* Beginning of gmo testing   */
    double *A;
    int m = 3; // number of rows of the matrix
    int n = 4; //number of columns of the matrix
    int numofCol;  //use to track how many columns the process has
    double finalAnser[m];
    double subValue; // variable to store sub value receive from each process
    // first, we need to allocate memory dynamically to store the matrix.
    A = (double *) malloc( 3 * 4 * sizeof(double));
    assert(A != NULL);
    genMatrix( 3, 4, A );
    printf(" the 2nd 1st number is: %f \n", *(A+4) );
    double *b;
    b = (double *) malloc( 4 * sizeof(double));
    assert(b!=NULL);
    genVector(4,b);
    double *AA, *bb, *result;
    int numofEntries;
    AA = (double *) malloc( 3 * 4 * sizeof(double));
    bb = (double *) malloc( 4 * sizeof(double));
    result = (double *) malloc( 45 * sizeof(double));
    assert( AA!=NULL && bb!=NULL);
    /* Start up MPI */
    MPI_Init(&argc, &argv);
    /* Find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    productDistributor( A, b, 3, 4, &numofEntries, AA, bb);
    internalProcessMultiplication( AA, bb, 3, numofEntries, result);
    
    // Now we need to send each of the result array to Master process so that to
    // product the final result
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank != 0 )
    {
        // each process will be sending out results based on row, this would ease the work for
        // receiver
        for( int i=0; i < m; i++ )
        {
            for (int j = 0; j < numofEntries; j++)
            {
                MPI_Send( (result + j*m + i), 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
            }
        }
        printf(" process is %i   gmo testing on numofEntries: %i \n", my_rank, numofEntries );
    }
    else  // this is master process
    {
        // first, we push the results in master process into final answer
        for (int i=0; i < m; i++)
        {
            for (int j=0; j < numofEntries; j++)
            {
                finalAnser[i] += *(result + j*m + i);
            }
        }
        // we need to receive each of message carefully here
        for (source = 1; source < p; source++)
        {
            // Since each of process could send out different number of messages depends on # of column it
            // contains, we have to figure out the number of entries. numofEntries is no valid here, it depends
            // on which process it is on.
            if (source < n%p)
            {
                numofCol = n/p + 1;
            }
            else
            {
                numofCol = n/p;
            }
            for (int i=0; i < m; i++)
            {
                for (int j=0; j < numofCol; j++)
                {
                    MPI_Recv(&subValue, 1, MPI_DOUBLE, source, tag,MPI_COMM_WORLD, &status);
                    printf("gmo received value is: %f \n", subValue);
                    finalAnser[i] += subValue;
                }
            }
        }
    }
    
    printf("final answer[0]: %f \n", finalAnser[0]);
    printf("final answer[1]: %f \n", finalAnser[1]);
    printf("final answer[2]: %f \n", finalAnser[2]);
    /* Shut down MPI */
    MPI_Finalize();
    
    return 0;
}
