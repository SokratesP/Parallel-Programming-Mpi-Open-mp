#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include<sys/wait.h>
#include<unistd.h>
#include <omp.h>

#define SIZE 16
#define UP    0
#define DOWN  1
#define LEFT  2
#define RIGHT 3
#define UPLEFT 4
#define UPRIGHT  5
#define DOWNLEFT  6
#define DOWNRIGHT 7

void main(int argc, char *argv[])  
{
  int numtasks ,i;
  int numOfPro;
  double start, finish, loc_elapsed, elapsed;
  int dimension=120;
  omp_set_num_threads(1);

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  numOfPro= (int)sqrt(numtasks);
  dimension= 1920/numtasks;
  int cont=0;
  for (i = 1; i < 12; ++i)
  {
    if(numtasks == i*i)
    {
      cont=1;
      break;
    }
  }

  if (cont == 1) 
  {
    int rank, source, dest, outbuf,j, tag=1, 
    inbuf[8]={MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL}, 
    periods[2]={0,0}, reorder=0, coords[2];
    int nbrs[8];
    int foundS=0;
    //time_t t;
    srand(time(NULL));
    //filtrer kai to normalization
    int filter[3][3]={{1,2,1},{2,4,2},{1,2,1}};
    int normalize=0;
    for (i = 0; i < 3; ++i)
    {
      normalize+=filter[i][0]+filter[i][1]+filter[i][2];
    }

    for (i = 0; i < 8; ++i)
    {
      nbrs[i] = -1;
    }
    MPI_Request reqs[16];
    MPI_Status stats[16];
    MPI_Comm cartcomm;  // required variable
    int dims[2]={numOfPro,numOfPro}; 
    // create cartesian virtual topology, get rank, coordinates, neighbor ranks
    //cartcomm = MPI_COMM_WORLD;


    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm);

    MPI_Barrier(cartcomm);
    start = MPI_Wtime();

    MPI_Comm_rank(cartcomm, &rank);
    MPI_Cart_coords(cartcomm, rank, 2, coords);

    MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UP], &nbrs[DOWN]);
    MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);

    //poia diergasia einai kai toytw diagvnious geitones
    outbuf = rank;
    if (rank==0)
    {
      nbrs[DOWNRIGHT]=nbrs[DOWN]+1;
    }else if (rank==numOfPro-1)
    {
      nbrs[DOWNLEFT]=nbrs[DOWN]-1;
    }else if (rank==numtasks - numOfPro)
    {
      nbrs[UPRIGHT]=nbrs[UP]+1;
    }else if (rank==numtasks-1)
    {
       nbrs[UPLEFT]=nbrs[UP]-1;
    }else if (rank>0 && rank < numOfPro-1)
    {
      nbrs[DOWNLEFT]=nbrs[DOWN]-1;
      nbrs[DOWNRIGHT]=nbrs[DOWN]+1;
    }else if ((rank > numtasks - numOfPro) && (rank < numtasks-1) )
    {
      nbrs[UPLEFT]=nbrs[UP]-1;
      nbrs[UPRIGHT]=nbrs[UP]+1;
    }else if (rank % numOfPro ==0)
    {
      //printf("Mphka gia rank %d\n",rank );
      nbrs[UPRIGHT]=nbrs[UP]+1;
      nbrs[DOWNRIGHT]=nbrs[DOWN]+1;
    }else if (rank % numOfPro == numOfPro -1 )
    {
      nbrs[UPLEFT]=nbrs[UP]-1;
      nbrs[DOWNLEFT]=nbrs[DOWN]-1;
    }else
    {
      nbrs[UPLEFT]=nbrs[UP]-1;
      nbrs[UPRIGHT]=nbrs[UP]+1;
      nbrs[DOWNLEFT]=nbrs[DOWN]-1;
      nbrs[DOWNRIGHT]=nbrs[DOWN]+1;
    }

    // exchange data (rank) with 8 neighbors
    for (i=0; i<8; i++) 
    {
      dest = nbrs[i];
      source = nbrs[i];
      MPI_Isend(&outbuf, 1, MPI_INT, dest, tag, MPI_COMM_WORLD, &reqs[i]);
      MPI_Irecv(&inbuf[i], 1, MPI_INT, source, tag, MPI_COMM_WORLD, &reqs[i+8]);
    }
    MPI_Waitall(16, reqs, stats);
   
 
    printf("rank= %d                  inbuf(u,d,l,r)= %d %d %d %d\n",
           rank,inbuf[UP],inbuf[DOWN],inbuf[LEFT],inbuf[RIGHT]);
    printf("rank= %d                  inbuf(UPLEFT,UPRIGHT,DOWNLEFT,DOWNRIGHT)= %d %d %d %d\n",
           rank,inbuf[UPLEFT],inbuf[UPRIGHT],inbuf[DOWNLEFT],inbuf[DOWNRIGHT]);         
    
    //periexei panti oti epejergazomaste
    int** matrix = (int **)malloc((dimension+2)*sizeof(int*));
    for (i = 0; i <dimension+2; ++i)
    {
      matrix[i]=(int *)malloc((dimension+2)*sizeof(int));
    }

    //periexei to apotelesma
    int** newMatrix = (int **)malloc((dimension)*sizeof(int*));
    for (i = 0; i <dimension; ++i)
    {
      newMatrix[i]=(int *)malloc((dimension)*sizeof(int));
    }

    //gemisma pinaka kai placeholder aspro gia tiw timew poy ua eruoyn
    for (i = 0; i < dimension+2; i++)
    {
      for (j = 0; j < dimension+2; j++)
      {

        if (i==0 || i==dimension+1 || j==0 || j == dimension+1 )
        {
          matrix[i][j]= 255;
        }
        else{
          matrix[i][j]=rand() % 256;
        }
      }
    }

   //------------------------------------------------------------------------------------------------//    
    //flags termatismouy
    int flag,prevFlag;

    //orizoyme ta noumera pou 8a steilei     
    int* upMatrix=(int *)malloc((dimension)*sizeof(int));
    int* downMatrix=(int *)malloc((dimension)*sizeof(int));
    int* leftMatrix=(int *)malloc((dimension)*sizeof(int));
    int* rightMatrix=(int *)malloc((dimension)*sizeof(int));
    //orizoume auta poy 8a laboume 

    int* recvupMatrix=(int *)malloc((dimension)*sizeof(int));
    int* recvdownMatrix=(int *)malloc((dimension)*sizeof(int));
    int* recvleftMatrix=(int *)malloc((dimension)*sizeof(int));
    int* recvrightMatrix=(int *)malloc((dimension)*sizeof(int));
    int s;

    for (s = 0; s < 1000; ++s)
    { 
      int upLeftNum=matrix[1][1];
      int upRightNum=matrix[1][dimension];
      int downLeftNum=matrix[dimension][1];
      int downRightNum=matrix[dimension][dimension];

      int recvupLeftNum=0;
      int recvupRightNum=0;
      int recvdownLeftNum=0;
      int recvdownRightNum=0;

      for (i = 1; i < dimension+1; ++i)
      {
        upMatrix[i-1]= matrix[1][i];
      }

      for (i = 1; i < dimension+1; ++i)
      {
        downMatrix[i-1]= matrix[dimension][i];
      }

      for (i = 1; i < dimension+1; ++i)
      {
        leftMatrix[i-1]= matrix[i][1];
      }

      for (i = 1; i < dimension+1; ++i)
      {
        rightMatrix[i-1]= matrix[i][dimension];
      }

      MPI_Request newReqs[16];
      MPI_Status newStats[16];
      //stelnei
      MPI_Isend(&upLeftNum, 1, MPI_INT, inbuf[UPLEFT], tag, MPI_COMM_WORLD, &newReqs[0]);
      MPI_Isend(&upMatrix[0], dimension, MPI_INT, inbuf[UP], tag, MPI_COMM_WORLD, &newReqs[1]);
      MPI_Isend(&upRightNum, 1, MPI_INT, inbuf[UPRIGHT], tag, MPI_COMM_WORLD, &newReqs[2]);
      MPI_Isend(&leftMatrix[0], dimension, MPI_INT, inbuf[LEFT], tag, MPI_COMM_WORLD, &newReqs[3]);
      MPI_Isend(&rightMatrix[0], dimension, MPI_INT, inbuf[RIGHT], tag, MPI_COMM_WORLD, &newReqs[4]);
      MPI_Isend(&downLeftNum, 1, MPI_INT, inbuf[DOWNLEFT], tag, MPI_COMM_WORLD, &newReqs[5]);
      MPI_Isend(&downMatrix[0], dimension, MPI_INT, inbuf[DOWN], tag, MPI_COMM_WORLD, &newReqs[6]);
      MPI_Isend(&downRightNum, 1, MPI_INT, inbuf[DOWNRIGHT], tag, MPI_COMM_WORLD, &newReqs[7]);
    //------------------------------------------------------------------------------------------------//

      //pairnei
      MPI_Irecv(&recvdownRightNum, 1, MPI_INT, inbuf[DOWNRIGHT], tag, MPI_COMM_WORLD, &newReqs[8]);   
      MPI_Irecv(&recvdownMatrix[0], dimension, MPI_INT, inbuf[DOWN], tag, MPI_COMM_WORLD, &newReqs[9]);
      MPI_Irecv(&recvdownLeftNum, 1, MPI_INT, inbuf[DOWNLEFT], tag, MPI_COMM_WORLD, &newReqs[10]);
      MPI_Irecv(&recvrightMatrix[0], dimension, MPI_INT, inbuf[RIGHT], tag, MPI_COMM_WORLD, &newReqs[11]);
      MPI_Irecv(&recvleftMatrix[0], dimension, MPI_INT, inbuf[LEFT], tag, MPI_COMM_WORLD, &newReqs[12]);
      MPI_Irecv(&recvupRightNum, 1, MPI_INT, inbuf[UPRIGHT], tag, MPI_COMM_WORLD, &newReqs[13]);
      MPI_Irecv(&recvupMatrix[0], dimension, MPI_INT, inbuf[UP], tag, MPI_COMM_WORLD, &newReqs[14]);
      MPI_Irecv(&recvupLeftNum, 1, MPI_INT, inbuf[UPLEFT], tag, MPI_COMM_WORLD, &newReqs[15]);

      //kanei to filtro sto eswteriko tou pnaka pou den xreiazetai kati apo ta recv
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      #pragma omp parallel for schedule(auto) 
      for (i = 2; i < dimension; i++)
      {
        for (j = 2; j < dimension; j++)
        {
          newMatrix[i-1][j-1]=((((matrix[i-1][j-1]*filter[0][0])+
                                  (matrix[i-1][j]*filter[0][1])+
                                  (matrix[i-1][j+1]*filter[0][2])+
                                  (matrix[i][j-1]*filter[1][0])+
                                  (matrix[i][j]*filter[1][1])+
                                  (matrix[i][j+1]*filter[1][2])+
                                  (matrix[i+1][j-1]*filter[2][0])+
                                  (matrix[i+1][j]*filter[2][1])+
                                  (matrix[i+1][j+1]*filter[2][2])))/normalize);
        }
      }
      //perimenei ta recv kai analoga me to poio ua labei, to bazei sthn uesh tou
      int recvIndex,k;
      //int pragmatemp;
      for (k = 0; k < 16; ++k)
      {
        MPI_Waitany(16, newReqs,&recvIndex,newStats);
        if (recvIndex == 8)
        {
          matrix[dimension+1][dimension+1] = recvdownRightNum;
        }
        else if (recvIndex == 9)
        {
          #pragma omp parallel for schedule(auto) 
          for (i = 1; i < dimension+1 ; ++i)
          {
            //int tid = omp_get_thread_num();
            //int total = omp_get_num_threads();

            //printf("Hello world from thread %d of %d\n", tid, total);
            matrix[dimension+1][i]=recvdownMatrix[i-1];
          }
          #pragma omp parallel for schedule(auto) 
          for (j = 2; j < dimension ; ++j)
          {
            if ((rank <= numtasks-1) && (rank >= (numtasks - numOfPro ) ) )
            {
              newMatrix[dimension-1][j-1]=((((matrix[dimension-1][j-1]*filter[0][0])+
                                  (matrix[dimension-1][j]*filter[0][1])+
                                  (matrix[dimension-1][j+1]*filter[0][2])+
                                  (matrix[dimension][j-1]*filter[1][0])+
                                  (matrix[dimension][j]*filter[1][1])+
                                  (matrix[dimension][j+1]*filter[1][2])+
                                  (matrix[dimension][j]*filter[2][0])+
                                  (matrix[dimension][j]*filter[2][1])+
                                  (matrix[dimension][j]*filter[2][2])))/normalize);
            }else
            {
              newMatrix[dimension-1][j-1]=((((matrix[dimension-1][j-1]*filter[0][0])+
                                  (matrix[dimension-1][j]*filter[0][1])+
                                  (matrix[dimension-1][j+1]*filter[0][2])+
                                  (matrix[dimension][j-1]*filter[1][0])+
                                  (matrix[dimension][j]*filter[1][1])+
                                  (matrix[dimension][j+1]*filter[1][2])+
                                  (matrix[dimension+1][j-1]*filter[2][0])+
                                  (matrix[dimension+1][j]*filter[2][1])+
                                  (matrix[dimension+1][j+1]*filter[2][2])))/normalize);
            }
            
          }
        }
        else if (recvIndex == 10)
        {
          matrix[dimension+1][0] = recvdownLeftNum;
        }
        else if (recvIndex == 11)
        {
          #pragma omp parallel for schedule(auto)
          for (i = 1; i < dimension+1 ; ++i)
          {
            matrix[i][dimension+1]=recvrightMatrix[i-1];
          }
          #pragma omp parallel for schedule(auto)
          for (i = 2; i < dimension ; ++i)
          {
            if (rank % numOfPro == numOfPro -1 )
            {
              newMatrix[i-1][dimension-1]=((((matrix[i-1][dimension-1]*filter[0][0])+
                                          (matrix[i-1][dimension]*filter[0][1])+
                                          (matrix[i][dimension]*filter[0][2])+
                                          (matrix[i][dimension-1]*filter[1][0])+
                                          (matrix[i][dimension]*filter[1][1])+
                                          (matrix[i][dimension]*filter[1][2])+
                                          (matrix[i+1][dimension-1]*filter[2][0])+
                                          (matrix[i+1][dimension]*filter[2][1])+
                                          (matrix[i][dimension]*filter[2][2])))/normalize);
            }else
            {
              newMatrix[i-1][dimension-1]=((((matrix[i-1][dimension-1]*filter[0][0])+
                                          (matrix[i-1][dimension]*filter[0][1])+
                                          (matrix[i-1][dimension+1]*filter[0][2])+
                                          (matrix[i][dimension-1]*filter[1][0])+
                                          (matrix[i][dimension]*filter[1][1])+
                                          (matrix[i][dimension+1]*filter[1][2])+
                                          (matrix[i+1][dimension-1]*filter[2][0])+
                                          (matrix[i+1][dimension]*filter[2][1])+
                                          (matrix[i+1][dimension+1]*filter[2][2])))/normalize);
            }
            
          }
        }
        else if (recvIndex == 12)
        {
          #pragma omp parallel for schedule(auto)
          for (i = 1; i < dimension+1 ; ++i)
          {
            matrix[i][0]=recvleftMatrix[i-1];
          }
          #pragma omp parallel for schedule(auto)
          for (i = 2; i < dimension ; ++i)
          {
            if (rank % numOfPro == 0 )
            {
              newMatrix[i-1][0]=((((matrix[i][1]*filter[0][0])+
                                  (matrix[i-1][1]*filter[0][1])+
                                  (matrix[i-1][2]*filter[0][2])+
                                  (matrix[i][1]*filter[1][0])+
                                  (matrix[i][1]*filter[1][1])+
                                  (matrix[i][2]*filter[1][2])+
                                  (matrix[i][1]*filter[2][0])+
                                  (matrix[i+1][1]*filter[2][1])+
                                  (matrix[i+1][2]*filter[2][2])))/normalize);
            }else
            {
              newMatrix[i-1][0]=((((matrix[i-1][0]*filter[0][0])+
                                  (matrix[i-1][1]*filter[0][1])+
                                  (matrix[i-1][2]*filter[0][2])+
                                  (matrix[i][0]*filter[1][0])+
                                  (matrix[i][1]*filter[1][1])+
                                  (matrix[i][2]*filter[1][2])+
                                  (matrix[i+1][0]*filter[2][0])+
                                  (matrix[i+1][1]*filter[2][1])+
                                  (matrix[i+1][2]*filter[2][2])))/normalize);
            }
            
          }
        }
        else if (recvIndex == 13)
        {
          matrix[0][dimension+1] = recvupRightNum;
        }
        else if (recvIndex == 14)
        {
          #pragma omp parallel for schedule(auto)
          for (i = 1; i < dimension+1 ; ++i)
          {
            matrix[0][i]=recvupMatrix[i-1];
          }
          #pragma omp parallel for schedule(auto)
          for (j = 2; j < dimension ; ++j)
          {
            if ((rank <= numOfPro-1) && (rank >= 0 ) )
            {
              newMatrix[0][j-1]=((((matrix[1][j]*filter[0][0])+
                                  (matrix[1][j]*filter[0][1])+
                                  (matrix[1][j]*filter[0][2])+
                                  (matrix[1][j-1]*filter[1][0])+
                                  (matrix[1][j]*filter[1][1])+
                                  (matrix[1][j+1]*filter[1][2])+
                                  (matrix[2][j-1]*filter[2][0])+
                                  (matrix[2][j]*filter[2][1])+
                                  (matrix[2][j+1]*filter[2][2])))/normalize);
            }else
            {
              newMatrix[0][j-1]=((((matrix[0][j-1]*filter[0][0])+
                                  (matrix[0][j]*filter[0][1])+
                                  (matrix[0][j+1]*filter[0][2])+
                                  (matrix[1][j-1]*filter[1][0])+
                                  (matrix[1][j]*filter[1][1])+
                                  (matrix[1][j+1]*filter[1][2])+
                                  (matrix[2][j-1]*filter[2][0])+
                                  (matrix[2][j]*filter[2][1])+
                                  (matrix[2][j+1]*filter[2][2])))/normalize);
            }
            
          }
        }
        else if (recvIndex == 15)
        {
          matrix[0][0]=recvupLeftNum;
        }
      }

      //diagwnia stoixeia afou exei parei ola ta stoixeia
      //(giati einai pio grhogro apo to na baleis 8 flag kai 16 elegxous)
      if (rank==0)
      {
        newMatrix[0][0]=((((matrix[1][1]*filter[0][0])+
                          (matrix[1][1]*filter[0][1])+
                          (matrix[1][1]*filter[0][2])+
                          (matrix[1][1]*filter[1][0])+
                          (matrix[1][1]*filter[1][1])+
                          (matrix[1][2]*filter[1][2])+
                          (matrix[1][1]*filter[2][0])+
                          (matrix[2][1]*filter[2][1])+
                          (matrix[2][2]*filter[2][2])))/normalize);

        newMatrix[0][dimension-1]=((((matrix[1][dimension]*filter[0][0])+
                                (matrix[1][dimension]*filter[0][1])+
                                (matrix[1][dimension]*filter[0][2])+
                                (matrix[1][dimension-1]*filter[1][0])+
                                (matrix[1][dimension]*filter[1][1])+
                                (matrix[1][dimension+1]*filter[1][2])+
                                (matrix[2][dimension-1]*filter[2][0])+
                                (matrix[2][dimension]*filter[2][1])+
                                (matrix[2][dimension+1]*filter[2][2])))/normalize);

        newMatrix[dimension-1][0]=((((matrix[dimension][1]*filter[0][0])+
                                  (matrix[dimension-1][1]*filter[0][1])+
                                  (matrix[dimension-1][2]*filter[0][2])+
                                  (matrix[dimension][1]*filter[1][0])+
                                  (matrix[dimension][1]*filter[1][1])+
                                  (matrix[dimension][2]*filter[1][2])+
                                  (matrix[dimension][1]*filter[2][0])+
                                  (matrix[dimension+1][1]*filter[2][1])+
                                  (matrix[dimension+1][2]*filter[2][2])))/normalize);

        newMatrix[dimension-1][dimension-1]=((((matrix[dimension-1][dimension-1]*filter[0][0])+
                                          (matrix[dimension-1][dimension]*filter[0][1])+
                                          (matrix[dimension-1][dimension+1]*filter[0][2])+
                                          (matrix[dimension][dimension-1]*filter[1][0])+
                                          (matrix[dimension][dimension]*filter[1][1])+
                                          (matrix[dimension][dimension+1]*filter[1][2])+
                                          (matrix[dimension+1][dimension-1]*filter[2][0])+
                                          (matrix[dimension+1][dimension]*filter[2][1])+
                                          (matrix[dimension+1][dimension+1]*filter[2][2])))/normalize);
      }
      else if (rank==numOfPro-1)
      {
        newMatrix[0][0]=((((matrix[1][1]*filter[0][0])+
                          (matrix[1][1]*filter[0][1])+
                          (matrix[1][1]*filter[0][2])+
                          (matrix[1][0]*filter[1][0])+
                          (matrix[1][1]*filter[1][1])+
                          (matrix[1][2]*filter[1][2])+
                          (matrix[2][0]*filter[2][0])+
                          (matrix[2][1]*filter[2][1])+
                          (matrix[2][2]*filter[2][2])))/normalize);

        newMatrix[0][dimension-1]=((((matrix[1][dimension]*filter[0][0])+
                                (matrix[1][dimension]*filter[0][1])+
                                (matrix[1][dimension]*filter[0][2])+
                                (matrix[1][dimension-1]*filter[1][0])+
                                (matrix[1][dimension]*filter[1][1])+
                                (matrix[1][dimension]*filter[1][2])+
                                (matrix[2][dimension-1]*filter[2][0])+
                                (matrix[2][dimension]*filter[2][1])+
                                (matrix[1][dimension]*filter[2][2])))/normalize);

        newMatrix[dimension-1][0]=((((matrix[dimension-1][0]*filter[0][0])+
                                  (matrix[dimension-1][1]*filter[0][1])+
                                  (matrix[dimension-1][2]*filter[0][2])+
                                  (matrix[dimension][0]*filter[1][0])+
                                  (matrix[dimension][1]*filter[1][1])+
                                  (matrix[dimension][2]*filter[1][2])+
                                  (matrix[dimension+1][0]*filter[2][0])+
                                  (matrix[dimension+1][1]*filter[2][1])+
                                  (matrix[dimension+1][2]*filter[2][2])))/normalize);

        newMatrix[dimension-1][dimension-1]=((((matrix[dimension-1][dimension-1]*filter[0][0])+
                                          (matrix[dimension-1][dimension]*filter[0][1])+
                                          (matrix[dimension][dimension]*filter[0][2])+
                                          (matrix[dimension][dimension-1]*filter[1][0])+
                                          (matrix[dimension][dimension]*filter[1][1])+
                                          (matrix[dimension][dimension]*filter[1][2])+
                                          (matrix[dimension+1][dimension-1]*filter[2][0])+
                                          (matrix[dimension+1][dimension]*filter[2][1])+
                                          (matrix[dimension][dimension]*filter[2][2])))/normalize);
      }
      else if (rank == numtasks - numOfPro)
      {
        newMatrix[0][0]=((((matrix[1][1]*filter[0][0])+
                          (matrix[0][1]*filter[0][1])+
                          (matrix[0][2]*filter[0][2])+
                          (matrix[1][1]*filter[1][0])+
                          (matrix[1][1]*filter[1][1])+
                          (matrix[1][2]*filter[1][2])+
                          (matrix[1][1]*filter[2][0])+
                          (matrix[2][1]*filter[2][1])+
                          (matrix[2][2]*filter[2][2])))/normalize);

        newMatrix[0][dimension-1]=((((matrix[0][dimension-1]*filter[0][0])+
                                (matrix[0][dimension]*filter[0][1])+
                                (matrix[0][dimension+1]*filter[0][2])+
                                (matrix[1][dimension-1]*filter[1][0])+
                                (matrix[1][dimension]*filter[1][1])+
                                (matrix[1][dimension+1]*filter[1][2])+
                                (matrix[2][dimension-1]*filter[2][0])+
                                (matrix[2][dimension]*filter[2][1])+
                                (matrix[2][dimension+1]*filter[2][2])))/normalize);

        newMatrix[dimension-1][0]=((((matrix[dimension][1]*filter[0][0])+
                                  (matrix[dimension-1][1]*filter[0][1])+
                                  (matrix[dimension-1][2]*filter[0][2])+
                                  (matrix[dimension][0]*filter[1][0])+
                                  (matrix[dimension][1]*filter[1][1])+
                                  (matrix[dimension][2]*filter[1][2])+
                                  (matrix[dimension][1]*filter[2][0])+
                                  (matrix[dimension][1]*filter[2][1])+
                                  (matrix[dimension][1]*filter[2][2])))/normalize);

        newMatrix[dimension-1][dimension-1]=((((matrix[dimension-1][dimension-1]*filter[0][0])+
                                          (matrix[dimension-1][dimension]*filter[0][1])+
                                          (matrix[dimension-1][dimension+1]*filter[0][2])+
                                          (matrix[dimension][dimension-1]*filter[1][0])+
                                          (matrix[dimension][dimension]*filter[1][1])+
                                          (matrix[dimension][dimension+1]*filter[1][2])+
                                          (matrix[dimension][dimension]*filter[2][0])+
                                          (matrix[dimension][dimension]*filter[2][1])+
                                          (matrix[dimension][dimension]*filter[2][2])))/normalize);
      }
      else if (rank == numtasks-1)
      {
        newMatrix[0][0]=((((matrix[0][0]*filter[0][0])+
                          (matrix[0][1]*filter[0][1])+
                          (matrix[0][2]*filter[0][2])+
                          (matrix[1][0]*filter[1][0])+
                          (matrix[1][1]*filter[1][1])+
                          (matrix[1][2]*filter[1][2])+
                          (matrix[2][0]*filter[2][0])+
                          (matrix[2][1]*filter[2][1])+
                          (matrix[2][2]*filter[2][2])))/normalize);

        newMatrix[0][dimension-1]=((((matrix[0][dimension-1]*filter[0][0])+
                                (matrix[0][dimension]*filter[0][1])+
                                (matrix[1][dimension]*filter[0][2])+
                                (matrix[1][dimension-1]*filter[1][0])+
                                (matrix[1][dimension]*filter[1][1])+
                                (matrix[1][dimension]*filter[1][2])+
                                (matrix[2][dimension-1]*filter[2][0])+
                                (matrix[2][dimension]*filter[2][1])+
                                (matrix[1][dimension]*filter[2][2])))/normalize);

        newMatrix[dimension-1][0]=((((matrix[dimension-1][0]*filter[0][0])+
                                  (matrix[dimension-1][1]*filter[0][1])+
                                  (matrix[dimension-1][2]*filter[0][2])+
                                  (matrix[dimension][0]*filter[1][0])+
                                  (matrix[dimension][1]*filter[1][1])+
                                  (matrix[dimension][2]*filter[1][2])+
                                  (matrix[dimension][1]*filter[2][0])+
                                  (matrix[dimension][1]*filter[2][1])+
                                  (matrix[dimension][1]*filter[2][2])))/normalize);

        newMatrix[dimension-1][dimension-1]=((((matrix[dimension-1][dimension-1]*filter[0][0])+
                                          (matrix[dimension-1][dimension]*filter[0][1])+
                                          (matrix[dimension][dimension]*filter[0][2])+
                                          (matrix[dimension][dimension-1]*filter[1][0])+
                                          (matrix[dimension][dimension]*filter[1][1])+
                                          (matrix[dimension][dimension]*filter[1][2])+
                                          (matrix[dimension][dimension]*filter[2][0])+
                                          (matrix[dimension][dimension]*filter[2][1])+
                                          (matrix[dimension][dimension]*filter[2][2])))/normalize);
      }else
      {
        newMatrix[0][0]=((((matrix[0][0]*filter[0][0])+
                          (matrix[0][1]*filter[0][1])+
                          (matrix[0][2]*filter[0][2])+
                          (matrix[1][0]*filter[1][0])+
                          (matrix[1][1]*filter[1][1])+
                          (matrix[1][2]*filter[1][2])+
                          (matrix[2][0]*filter[2][0])+
                          (matrix[2][1]*filter[2][1])+
                          (matrix[2][2]*filter[2][2])))/normalize);

        newMatrix[0][dimension-1]=((((matrix[0][dimension-1]*filter[0][0])+
                                (matrix[0][dimension]*filter[0][1])+
                                (matrix[0][dimension+1]*filter[0][2])+
                                (matrix[1][dimension-1]*filter[1][0])+
                                (matrix[1][dimension]*filter[1][1])+
                                (matrix[1][dimension+1]*filter[1][2])+
                                (matrix[2][dimension-1]*filter[2][0])+
                                (matrix[2][dimension]*filter[2][1])+
                                (matrix[2][dimension+1]*filter[2][2])))/normalize);

        newMatrix[dimension-1][0]=((((matrix[dimension-1][0]*filter[0][0])+
                                  (matrix[dimension-1][1]*filter[0][1])+
                                  (matrix[dimension-1][2]*filter[0][2])+
                                  (matrix[dimension][0]*filter[1][0])+
                                  (matrix[dimension][1]*filter[1][1])+
                                  (matrix[dimension][2]*filter[1][2])+
                                  (matrix[dimension+1][0]*filter[2][0])+
                                  (matrix[dimension+1][1]*filter[2][1])+
                                  (matrix[dimension+1][2]*filter[2][2])))/normalize);

        newMatrix[dimension-1][dimension-1]=((((matrix[dimension-1][dimension-1]*filter[0][0])+
                                          (matrix[dimension-1][dimension]*filter[0][1])+
                                          (matrix[dimension-1][dimension+1]*filter[0][2])+
                                          (matrix[dimension][dimension-1]*filter[1][0])+
                                          (matrix[dimension][dimension]*filter[1][1])+
                                          (matrix[dimension][dimension+1]*filter[1][2])+
                                          (matrix[dimension+1][dimension-1]*filter[2][0])+
                                          (matrix[dimension+1][dimension]*filter[2][1])+
                                          (matrix[dimension+1][dimension+1]*filter[2][2])))/normalize);
      }
      
      int sum=0;
      //#pragma omp parallel for
      for (i = 0; i < dimension+2; i++)
      {
        for (j = 0; j < dimension+2; j++)
        {

          if (i==0 || i==dimension+1 || j==0 || j == dimension+1 )
          {
            matrix[i][j]= 255;
          }
          else{
            if(newMatrix[i-1][j-1]>=255)
            {
              matrix[i][j]=255;
            }else if (newMatrix[i-1][j-1]<=0)
            {
              //matrix[i][j]=0;
            }
            else
            {
              matrix[i][j]=newMatrix[i-1][j-1];
            }
            //#pragma omp critical  
            sum+=matrix[i][j];
          }
        }
      }
      flag=sum/(dimension*dimension);
      if (flag < prevFlag )
      {
        foundS = s;
      }
      prevFlag = flag;
    }
    //------------------------------------------------------------------------------------------------/

    //printf("Rank :::: %d s= %d     foundS=%d \n",rank ,s ,foundS);

    /*for (i = 0; i < dimension; ++i)
    {
      for (j = 0; j < dimension; ++j)
      {
        printf("%d  ", newMatrix[i][j]);
      }
      printf("\n");
    }*/
    finish = MPI_Wtime();
    loc_elapsed = finish-start;
    MPI_Reduce(&loc_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, cartcomm);

    if (rank == 0)
    {
      printf("Elapsed time = %e\n", elapsed);
    }
    free(upMatrix);
    free(downMatrix);
    free(rightMatrix);
    free(leftMatrix);
    free(recvupMatrix);
    free(recvdownMatrix);
    free(recvleftMatrix);
    free(recvrightMatrix);
    for (i = 0; i < dimension+2; ++i)
    {
      free(matrix[i]);
      if (i<dimension)
      {
        free(newMatrix[i]);
      }
    }
    free(matrix);
    free(newMatrix);
    MPI_Finalize();

  }
  else{
    printf("Number of tasks can not be divided. Terminating. :)\n");
    //MPI_Finalize();
   }
   
  
  exit(0);
}
