#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mpi.h"
#include <assert.h>
#define INTERVAL_FACTOR 5000 //This means we will devide the range into 50 sub ranges
#define TERMINATE_SIGNAL 1 //To signal work is done

unsigned long int findNextPrime( unsigned long int curNum );
unsigned long int findLargestPrimeGap( unsigned long int beginVal, unsigned long int endVal);
unsigned long int findLargestPrimeGapLast( unsigned long int beginVal, unsigned long int endVal);
void * Malloc(int num, int size, MPI_Comm comm, int rank)
{
        void *b = malloc(num*size);
        if( b== NULL)
        {
                fprintf(stderr," **** Process %d: ALLOC COULD NOT ALLOCATE %d Elements of size %d bytes \n", rank, num, size);
                MPI_Abort(comm, 1);
        }
        return b;
}

double elapsed_time;

int main(int argc, char* argv[])
{

       int      my_rank;
       int      tag = 1;       /* tag for messages     */
       int      ftag = 2; /*tag for final result comminication */
       int      numProcess;
       unsigned long int beginValue, endValue; // variables that stored the end values of the range
       char *initialValue, *lastValue;
       initialValue = argv[1];
       lastValue = argv[2];
       beginValue = strtoull (initialValue, &initialValue, 10); //convert str to llu
       endValue = strtoull( lastValue, &lastValue, 10);

       MPI_Init(&argc,&argv);
       MPI_Comm_size(MPI_COMM_WORLD, &numProcess);
       MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
       elapsed_time = - MPI_Wtime(); // for checking performance purpose

    if (numProcess != 1)
    {
       // code begining, assign process 0 to be the main process
       if( my_rank == 0 )
       {
           //Main process is responsible to assign work for other process to
           //find the largest gap within each range, and it also is responsible to
           //find gaps in a relatively smaller range, and find out the largest gap.
           MPI_Request reqs[numProcess-1];  // for communication check purposes
           MPI_Request reqsDone[numProcess-1];  // use to check if it is done!
           MPI_Status  status[numProcess-1];        /* status for receive   */
           unsigned long int eachProRange = ((endValue - beginValue)/INTERVAL_FACTOR)/(numProcess);
           int slaveProcessWork = numProcess/2;   //each slave is responsible for more work
           int mainProcessWork = 1;
           int maxGap = 0; // maximum prime gap for all ranges
           int curCount = 0; // counter used by main process
           // allocate an array to keep track on each assigned values
           int i;
           int flags[numProcess -1];   //define flags for communication purpose
           for (i = 0; i < (numProcess -1); i++)
           {
               flags[i] =1;
           }


           int totalPartitions = INTERVAL_FACTOR * (numProcess);  //number of total partitions divided
           //an array that holds the end bounds of each range
           unsigned long int * rangeArray = (unsigned long int * ) Malloc( totalPartitions+1, sizeof(unsigned long int ), MPI_COMM_WORLD, my_rank );
           int revBuffer[(numProcess -1)]; //buffer that main process use to collect max gap from slave
           unsigned long int sendBuffer[2]; //contains the boundaries of a range
           int signalBuffer[(numProcess -1)]; //buffer that main process use to collect signals

           // S0: divide all work into subranges, and stores them in rangePartition
           rangeArray[0] = beginValue;
           int endBound = (numProcess-1)* slaveProcessWork + 1;
           for( i = 1; i < endBound; i++)
           {
               rangeArray[i] = rangeArray[i-1] + eachProRange;
           }

           // first task distribution to other process, so that other process can start work
           // meanwhile we keep generating the other entries.
           for( i = 1; i < numProcess; i++)
           {
               sendBuffer[0] = rangeArray[curCount];
               sendBuffer[1] = rangeArray[curCount + slaveProcessWork];
               MPI_Isend(&sendBuffer, 2, MPI_UNSIGNED_LONG, i, tag, MPI_COMM_WORLD, &reqs[i] );
               curCount += slaveProcessWork;  //use to keep track on how many send requests had been sent
           }

           for( i = endBound; i < totalPartitions+1; i++)
           {
               rangeArray[i] = rangeArray[i-1] + eachProRange;
           }


           // we should change the last range to endValue to cover any missing range due to division
           rangeArray[totalPartitions] = endValue;
           // from now on, we should keep track on a while loop
           int Done = 0;
           int mainGap = 0;
           while( !Done)
           {
                // check if a send has been done based on the flag value
                for( i = 0; i< (numProcess -1 ); i++)
                {
                    if( flags[i] )
                    {
                        MPI_Irecv( &signalBuffer[i], 1, MPI_INT, i+1, tag, MPI_COMM_WORLD, &reqs[i]);
                    }
                }

                // check if we have received from slave process
                for( i = 0; i< numProcess-1; i++)
                {
                    MPI_Test(&reqs[i], &flags[i], &status[i]);
                }

                for( i = 0; i < (numProcess - 1); i++ )
                {
                    if( flags[i] )
                    {
                        // if we have received stuff, we also proceed to store the values, in the meantime, we can
                        // again send more stuff to the process to do work
                        if (curCount + slaveProcessWork > totalPartitions)
                        {
                            if (curCount != totalPartitions)
                            {
                                sendBuffer[0] = rangeArray[curCount];
                                sendBuffer[1] = rangeArray[totalPartitions]; //move the end bound to the final position
                                MPI_Isend(&sendBuffer, 2, MPI_UNSIGNED_LONG, i+1, tag, MPI_COMM_WORLD, &reqs[i] );
                                curCount = totalPartitions; //end of the array
                            }

                            //also, we should now send out signal to all slaves to terminate all operations
                            Done = 1;
                            unsigned long int endBuff[2];
                            endBuff[1] = TERMINATE_SIGNAL;
                            // send out message to all process to inform it is done
                            for( i = 1; i < numProcess; i++)
                            {
                                MPI_Isend( &endBuff, 2, MPI_UNSIGNED_LONG, i, tag, MPI_COMM_WORLD, &reqs[i-1] );
                                MPI_Irecv( &revBuffer[i-1], 1, MPI_INT, i, ftag, MPI_COMM_WORLD, &reqsDone[i-1]);
                                // MPI_Irecv( &revBuffer[i-1], 1, MPI_INT, i, ftag, MPI_COMM_WORLD, &reqsDone[i-1]);
                            }
                            break;
                        }
                        else
                        {
                            // normal routine of assigning work to slaves
                            sendBuffer[0] = rangeArray[curCount];
                            sendBuffer[1] = rangeArray[curCount + slaveProcessWork];
                            MPI_Isend(&sendBuffer, 2, MPI_UNSIGNED_LONG, i+1, tag, MPI_COMM_WORLD, &reqs[i] );
                            curCount += slaveProcessWork;
                        }
                    }
                }

                // we can make use of main process for computation
               if (curCount < totalPartitions)
               {
                   if ( (curCount + mainProcessWork) != totalPartitions)
                   {
                       mainGap = findLargestPrimeGap( rangeArray[curCount], rangeArray[curCount + mainProcessWork] );
                   }
                   else
                   {
                       mainGap = findLargestPrimeGapLast( rangeArray[curCount], rangeArray[curCount + mainProcessWork] );
                   }
                   curCount+= mainProcessWork;
                   if (mainGap > maxGap) maxGap = mainGap;
               }
            }  //end of while loop.

            // now we need to poll to receive all the largest gaps from slaves
           MPI_Waitall( numProcess-1, reqsDone, status);
           for (i=0; i< numProcess -1; i++)
           {
               if (revBuffer[i] > maxGap)
               {
                   maxGap = revBuffer[i];
               }
           }

           free(rangeArray);
           printf("THE MAXMUM GAP is %d \n", maxGap);
       }
       else  // other process
       {
                // S1: listen from main process for task
                unsigned long int subInterval[2];
                int flag=0;
                MPI_Status stat;
                MPI_Request req;
                int done = 0;
                MPI_Irecv( &subInterval, 2, MPI_UNSIGNED_LONG, 0, tag, MPI_COMM_WORLD, &req);
                int localMaxGap = 0;
                int processMaxGap = 0;
                int callCount =0; //gmo
                while ( !done)
                {
                    MPI_Test(&req, &flag, &stat);
                    // S2: compute task
                    if(flag && (subInterval[1] != TERMINATE_SIGNAL))
                    {
                        MPI_Isend( &localMaxGap, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &req ); //to signal main process for more work!
                        if (subInterval[1] != endValue)
                        {
                            localMaxGap = findLargestPrimeGap( subInterval[0], subInterval[1] );
                        }
                        else
                        {
                            localMaxGap = findLargestPrimeGapLast( subInterval[0], subInterval[1] );
                        }

                        if (localMaxGap > processMaxGap) processMaxGap = localMaxGap;

                        callCount++;
                        MPI_Irecv( subInterval, 2, MPI_UNSIGNED_LONG, 0, tag, MPI_COMM_WORLD, &req);
                    }
                    if ( subInterval[1] == TERMINATE_SIGNAL)
                    {
                        done = 1;
                    }
                }
                MPI_Isend( &processMaxGap, 1, MPI_INT, 0, ftag, MPI_COMM_WORLD, &req );
                printf("the max G at P: %d is %d, called %d times \n", my_rank, processMaxGap, callCount);
       }

        elapsed_time += MPI_Wtime();
        printf(" the time is: %f s , and P: %d \n", elapsed_time, my_rank);
    }
    else //only have one process
    {
        int maxGap = 0;
        maxGap = findLargestPrimeGapLast( beginValue, endValue );
        elapsed_time += MPI_Wtime();
        printf(" the time is: %f s , and P: %d \n", elapsed_time, my_rank);
        printf("THE MAXMUM GAP is %d \n", maxGap);
    }
       MPI_Finalize();
       return 0;
}

// This is just a wrapper function to find next prime number
unsigned long int findNextPrime( unsigned long int curNum )
{
        mpz_t nextPrime, curInput;
        mpz_init( nextPrime );
        mpz_init( curInput );
        unsigned long int result;
        // convert int to string
        char str[64];
        sprintf( str, "%lu", curNum);
        mpz_set_str ( curInput, str, 10 );
        mpz_nextprime( nextPrime, curInput );
        mpz_get_str ( str, 10, nextPrime );
        mpz_clear( nextPrime );
        mpz_clear( curInput );
	result = strtoull( str, &str, 10);
        return result;
}

/*
 * Function: findLargestPrimeGap( int beginVal, int endVal )
 * Description: This function finds the largest gap within an interval between
 *              beginVal and endVal, and it also sends the first and last
 *              prime numbers in this range.
 */
unsigned long int findLargestPrimeGap( unsigned long int beginVal, unsigned long int endVal)
{
        int maxGap = 0;
        int diff;
        unsigned long int curPrime, nextPrime;
        curPrime = findNextPrime( beginVal - 1 );

        //keep generating prime number to the array till hit the end val, please note that this also
        //includes the first prime number that pass the end limit, so that we dont need to worry about
        //prime gap between different interval.
        while( ( nextPrime = findNextPrime( curPrime ) ) <= endVal )
        {
            diff = nextPrime - curPrime;
            if( maxGap < diff ) maxGap = diff;
            curPrime = nextPrime;
        }
        diff = nextPrime - curPrime;
        if( maxGap < diff ) maxGap = diff;
        return maxGap;
}

/*
 * Function: findLargestPrimeGapLast( int beginVal, int endVal )
 * Description: Duplication of findlargestPrimeGap as C does not support default function
 *              This is used for last sub range
 */
unsigned long int findLargestPrimeGapLast( unsigned long int beginVal, unsigned long int endVal)
{
    int maxGap = 0;
    int diff;
    unsigned long int curPrime, nextPrime;
    curPrime = findNextPrime( beginVal - 1 );
    //keep generating prime number to the array till hit the end val
    while( ( nextPrime = findNextPrime( curPrime ) ) <= endVal )
    {
        diff = nextPrime - curPrime;
        if( maxGap < diff ) maxGap = diff;
        curPrime = nextPrime;
    }
    return maxGap;
}

