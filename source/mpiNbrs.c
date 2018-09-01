#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include<sys/wait.h>
#include<unistd.h>

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
  int numtasks, rank, source, dest, outbuf, i,j, tag=1, 
  inbuf[8]={MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL}, 
  dims[2]={4,4}, 
  periods[2]={0,0}, reorder=0, coords[2];
  int nbrs[8];
  //time_t t;
  srand(time(NULL));

  for (i = 0; i < 8; ++i)
  {
    nbrs[i] = -1;
  }

  int numOfPro;
  MPI_Request reqs[16];
  MPI_Status stats[16];
  MPI_Comm cartcomm;  // required variable

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  numOfPro= (int)sqrt(numtasks);

  //printf("%d\n",numOfPro );
  if (numtasks == SIZE) 
  {
      // create cartesian virtual topology, get rank, coordinates, neighbor ranks
      MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm);
      MPI_Comm_rank(cartcomm, &rank);
      MPI_Cart_coords(cartcomm, rank, 2, coords);

      MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UP], &nbrs[DOWN]);
      MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);

      //MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UPLEFT], &nbrs[DOWNLEFT]);
      //MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);
      printf("#A: rank= %d coords= %d %d  neighbors(u,d,l,r)= %d %d %d %d\n",
             rank,coords[0],coords[1],nbrs[UP],nbrs[DOWN],nbrs[LEFT],
             nbrs[RIGHT]);

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

      printf("rank= %d coords= %d %d  neighbors(UPLEFT,UPRIGHT,DOWNLEFT,DOWNRIGHT)= %d %d %d %d\n",
             rank,coords[0],coords[1],nbrs[UPLEFT],nbrs[UPRIGHT],nbrs[DOWNLEFT],
             nbrs[DOWNRIGHT]);

      // exchange data (rank) with 4 neighbors
      for (i=0; i<8; i++) {
         dest = nbrs[i];
         source = nbrs[i];
         MPI_Isend(&outbuf, 1, MPI_INT, dest, tag, 
                   MPI_COMM_WORLD, &reqs[i]);
         MPI_Irecv(&inbuf[i], 1, MPI_INT, source, tag, 
                   MPI_COMM_WORLD, &reqs[i+8]);
         }

      MPI_Waitall(16, reqs, stats);
   
      printf("rank= %d                  inbuf(u,d,l,r)= %d %d %d %d\n",
             rank,inbuf[UP],inbuf[DOWN],inbuf[LEFT],inbuf[RIGHT]);
      printf("rank= %d                  inbuf(UPLEFT,UPRIGHT,DOWNLEFT,DOWNRIGHT)= %d %d %d %d\n",
             rank,inbuf[UPLEFT],inbuf[UPRIGHT],inbuf[DOWNLEFT],inbuf[DOWNRIGHT]);         
      
      int** matrix = malloc((numOfPro+2)*sizeof(int*));
      for (i = 0; i <numOfPro+2; ++i)
      {
        matrix[i]=malloc((numOfPro+2)*sizeof(int));
      }
      printf("###AAAAAAAA######\n" );
      for (i = 0; i < numOfPro+2; i++)
      {
        for (j = 0; j < numOfPro+2; j++)
        {
          if (i==0 || i==numOfPro+1 || j==0 || j == numOfPro+1 )
          {
            matrix[i][j]= -4;
          }
          else{
            matrix[i][j]=rand() % 255;  
          }
        }
      }
      printf("###bbbbbbbbbbbb######\n" );
      sleep(2);
      printf("Rank : %d\n",rank );
      for (i = 0; i < numOfPro+2; ++i)
      {
        for (j = 0; j < numOfPro+2; ++j)
        {
          printf("%d  ", matrix[i][j]);
        }
        printf("\n");
      }


    }
   else{
      printf("Must specify %d processors. Terminating.\n",SIZE);
   }
   
  MPI_Finalize();
}
