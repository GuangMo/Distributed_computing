#include <stdio.h>
#include <string.h> 
#include <unistd.h>
#include "mpi.h" 
#include "genmatvec.c" 
#include "matvecres.c"

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
    int i, j;  
    for ( i = 0; i < n; i++)
    {
        if( my_rank == i%numP ) // the assgined process which takes on the corresponding ith column and ith entry
        {
            // we need to push m entries to array AA to represent the entire column of A, or we can convert
            // matrix A as coloum based matrix by how the way matrix was generated, this could improve performance significantly as it reduces lots of Cache miss.
            for ( j = 0; j < m; j++)
            {
                *(AA + (*numofEntries)*m + j) = *(A + i + (j*n));
            }
            *(bb + *numofEntries) = *(b + i);
            *numofEntries += 1;    // to signal this process now contains one more entry
        }
    }
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
    int i,j;
    for ( i=0; i < numofEntries; i++)
    {
        // multiplication for each column
        for ( j = 0; j < m; j++)
        {
            *(result + i*m + j) = (*(AA + i*m + j)) * (*(bb + i));
        }
    }
}

double elapsed_time; 

int main(int argc, char* argv[])
{
    int         my_rank;       /* rank of process      */
    int         p;             /* number of processes  */
    int         source;        /* rank of sender       */
    int         tag = 0;       /* tag for messages     */
    MPI_Status  status;        /* status for receive   */
    double *A, *b; // variables uses for global matrix and vector
    int m = atoi(argv[1]); // number of rows of the matrix
    int n = atoi(argv[2]); //number of columns of the matrix
    int numofCol;  //use to track how many columns the process has
    double *finalAnser; //array that uses to store the final answers after computation.
    finalAnser = (double *) malloc( m * sizeof(double));
    // first, we need to allocate memory dynamically to store the matrix.
    A = (double *) malloc( m * n * sizeof(double));
    assert(A != NULL);
    b = (double *) malloc( n * sizeof(double));
    assert(b!=NULL);
    genMatrix( m, n, A );       //generate matrix entries randomly
    genVector(n,b);             //generate vectore entries randomly
    double *AA, *bb, *result;   //arrays that stores local values for each process
    int numofEntries;           //variable that keep tracks on how many entries in each P
    /* Start up MPI */
    MPI_Init(&argc, &argv);
    /* Find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    //allocate memory for arrays to each process, which will be used to contain part of the
    //global matrix and vector, and the size will use maximum possible size
    int maxPossibleEntries = n/p + 1;
    AA = (double *) malloc( m * maxPossibleEntries * sizeof(double));
    bb = (double *) malloc( maxPossibleEntries * sizeof(double));
    result = (double *) malloc( m * maxPossibleEntries * sizeof(double));
    assert( AA!=NULL && bb!=NULL && result!=NULL);
    elapsed_time = - MPI_Wtime(); // for checking performance purpose
    productDistributor( A, b, m, n, &numofEntries, AA, bb); //Distribute the matrix to each P
    internalProcessMultiplication( AA, bb, m, numofEntries, result); //Computes the product for each P
    // Now we need to send each of the result array to Master process so that to
    // product the final result
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank != 0 )
    {
        int numofSend = m * numofEntries; //send all entries in the result array from each P
        MPI_Send( result, numofSend, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }
    else  // this is master process
    {
        // first, we push the results in master process into final answer
        int i, j;
        for ( i=0; i < m; i++)
        {
            for ( j=0; j < numofEntries; j++)
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
            
            double *subResults;
            subResults = (double *) malloc(m*numofEntries *sizeof(double)); //process zero will always holds largest amount of data, so use its values to allocate an array would be able to store all values from other process.
            assert(subResults != NULL); 
            MPI_Recv(subResults, numofCol*m, MPI_DOUBLE, source, tag,MPI_COMM_WORLD, &status);
            // now we need to decode the data from receive array
            int i, j;
            for(i=0; i<m; i++)
            {
                for ( j = 0; j < numofCol; j++)
                {
                    finalAnser[i] += subResults[i + j*m];
                }
            }
            free( subResults);
        }
        elapsed_time += MPI_Wtime();
        printf(" the time is: %f s \n", elapsed_time);
    }
    getResult( m, n, A, b, finalAnser);//let TA check my work
    free(A);free(b);free(AA);free(bb);
    free(result);free(finalAnser);
    /* Shut down MPI */
    MPI_Finalize();
    return 0;
}
