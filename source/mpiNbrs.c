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
  int foundS=0;
  //time_t t;
  srand(time(NULL));
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
    int** matrix = malloc((numOfPro+2)*sizeof(int*));
    for (i = 0; i <numOfPro+2; ++i)
    {
      matrix[i]=malloc((numOfPro+2)*sizeof(int));
    }

    //periexei to apotelesma
    int** newMatrix = malloc((numOfPro)*sizeof(int*));
    for (i = 0; i <numOfPro; ++i)
    {
      newMatrix[i]=malloc((numOfPro)*sizeof(int));
    }

    for (i = 0; i < numOfPro+2; i++)
    {
      for (j = 0; j < numOfPro+2; j++)
      {

        if (i==0 || i==numOfPro+1 || j==0 || j == numOfPro+1 )
        {
          matrix[i][j]= 255;
        }
        else{
          matrix[i][j]=rand() % 256;  
          //newMatrix[i-1][j-1] = 0;
        }
      }
    }

//------------------------------------------------------------------------------------------------//    
    int s;
    for (s = 0; s < 100; ++s)
    {      
      int upLeftNum=matrix[1][1];
      int upRightNum=matrix[1][numOfPro];
      int downLeftNum=matrix[numOfPro][1];
      int downRightNum=matrix[numOfPro][numOfPro];

      int* upMatrix=malloc((numOfPro)*sizeof(int));
      for (i = 1; i < numOfPro+1; ++i)
      {
        upMatrix[i-1]= matrix[1][i];
      }

      int* downMatrix=malloc((numOfPro)*sizeof(int));
      for (i = 1; i < numOfPro+1; ++i)
      {
        downMatrix[i-1]= matrix[numOfPro][i];
      }

      int* leftMatrix=malloc((numOfPro)*sizeof(int));
      for (i = 1; i < numOfPro+1; ++i)
      {
        leftMatrix[i-1]= matrix[i][1];
      }

      int* rightMatrix=malloc((numOfPro)*sizeof(int));
      for (i = 1; i < numOfPro+1; ++i)
      {
        rightMatrix[i-1]= matrix[i][numOfPro];
      }

      MPI_Request newReqs[16];
      MPI_Status newStats[16];

      MPI_Isend(&upLeftNum, 1, MPI_INT, inbuf[UPLEFT], tag, MPI_COMM_WORLD, &newReqs[0]);
      MPI_Isend(&upMatrix[0], numOfPro, MPI_INT, inbuf[UP], tag, MPI_COMM_WORLD, &newReqs[1]);
      MPI_Isend(&upRightNum, 1, MPI_INT, inbuf[UPRIGHT], tag, MPI_COMM_WORLD, &newReqs[2]);
      MPI_Isend(&leftMatrix[0], numOfPro, MPI_INT, inbuf[LEFT], tag, MPI_COMM_WORLD, &newReqs[3]);
      MPI_Isend(&rightMatrix[0], numOfPro, MPI_INT, inbuf[RIGHT], tag, MPI_COMM_WORLD, &newReqs[4]);
      MPI_Isend(&downLeftNum, 1, MPI_INT, inbuf[DOWNLEFT], tag, MPI_COMM_WORLD, &newReqs[5]);
      MPI_Isend(&downMatrix[0], numOfPro, MPI_INT, inbuf[DOWN], tag, MPI_COMM_WORLD, &newReqs[6]);
      MPI_Isend(&downRightNum, 1, MPI_INT, inbuf[DOWNRIGHT], tag, MPI_COMM_WORLD, &newReqs[7]);
//------------------------------------------------------------------------------------------------//
  
      int recvupLeftNum=0;
      int recvupRightNum=0;
      int recvdownLeftNum=0;
      int recvdownRightNum=0;
      int* recvupMatrix=malloc((numOfPro)*sizeof(int));
      int* recvdownMatrix=malloc((numOfPro)*sizeof(int));
      int* recvleftMatrix=malloc((numOfPro)*sizeof(int));
      int* recvrightMatrix=malloc((numOfPro)*sizeof(int));



      MPI_Irecv(&recvdownRightNum, 1, MPI_INT, inbuf[DOWNRIGHT], tag, MPI_COMM_WORLD, &newReqs[8]);   
      MPI_Irecv(&recvdownMatrix[0], numOfPro, MPI_INT, inbuf[DOWN], tag, MPI_COMM_WORLD, &newReqs[9]);
      MPI_Irecv(&recvdownLeftNum, 1, MPI_INT, inbuf[DOWNLEFT], tag, MPI_COMM_WORLD, &newReqs[10]);
      MPI_Irecv(&recvrightMatrix[0], numOfPro, MPI_INT, inbuf[RIGHT], tag, MPI_COMM_WORLD, &newReqs[11]);
      MPI_Irecv(&recvleftMatrix[0], numOfPro, MPI_INT, inbuf[LEFT], tag, MPI_COMM_WORLD, &newReqs[12]);
      MPI_Irecv(&recvupRightNum, 1, MPI_INT, inbuf[UPRIGHT], tag, MPI_COMM_WORLD, &newReqs[13]);
      MPI_Irecv(&recvupMatrix[0], numOfPro, MPI_INT, inbuf[UP], tag, MPI_COMM_WORLD, &newReqs[14]);
      MPI_Irecv(&recvupLeftNum, 1, MPI_INT, inbuf[UPLEFT], tag, MPI_COMM_WORLD, &newReqs[15]);

      for (i = 2; i < numOfPro; i++)
      {

        for (j = 2; j < numOfPro; j++)
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
      int recvIndex,k;
      for (k = 0; k < 16; ++k)
      {
        MPI_Waitany(16, newReqs,&recvIndex,newStats);
        if (recvIndex == 8)
        {
          matrix[numOfPro+1][numOfPro+1] = recvdownRightNum;
          //printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$mphka\n");
        }
        else if (recvIndex == 9)
        {
          for (i = 1; i < numOfPro+1 ; ++i)
          {
            matrix[numOfPro+1][i]=recvdownMatrix[i-1];
          }

          for (j = 2; j < numOfPro ; ++j)
          {
            newMatrix[numOfPro-1][j-1]=((((matrix[numOfPro-1][j-1]*filter[0][0])+
                                  (matrix[numOfPro-1][j]*filter[0][1])+
                                  (matrix[numOfPro-1][j+1]*filter[0][2])+
                                  (matrix[numOfPro][j-1]*filter[1][0])+
                                  (matrix[numOfPro][j]*filter[1][1])+
                                  (matrix[numOfPro][j+1]*filter[1][2])+
                                  (matrix[numOfPro+1][j-1]*filter[2][0])+
                                  (matrix[numOfPro+1][j]*filter[2][1])+
                                  (matrix[numOfPro+1][j+1]*filter[2][2])))/normalize);
          }
        }
        else if (recvIndex == 10)
        {
          matrix[numOfPro+1][0] = recvdownLeftNum;
        }
        else if (recvIndex == 11)
        {
          for (i = 1; i < numOfPro+1 ; ++i)
          {
            matrix[i][numOfPro+1]=recvrightMatrix[i-1];
          }

          for (i = 2; i < numOfPro ; ++i)
          {
            newMatrix[i-1][numOfPro-1]=((((matrix[i-1][numOfPro-1]*filter[0][0])+
                                          (matrix[i-1][numOfPro]*filter[0][1])+
                                          (matrix[i-1][numOfPro+1]*filter[0][2])+
                                          (matrix[i][numOfPro-1]*filter[1][0])+
                                          (matrix[i][numOfPro]*filter[1][1])+
                                          (matrix[i][numOfPro+1]*filter[1][2])+
                                          (matrix[i+1][numOfPro-1]*filter[2][0])+
                                          (matrix[i+1][numOfPro]*filter[2][1])+
                                          (matrix[i+1][numOfPro+1]*filter[2][2])))/normalize);
          }
        }
        else if (recvIndex == 12)
        {
          for (i = 1; i < numOfPro+1 ; ++i)
          {
            matrix[i][0]=recvleftMatrix[i-1];
          }

          for (i = 2; i < numOfPro ; ++i)
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
        else if (recvIndex == 13)
        {
          matrix[0][numOfPro+1] = recvupRightNum;
        }
        else if (recvIndex == 14)
        {
          for (i = 1; i < numOfPro+1 ; ++i)
          {
            matrix[0][i]=recvupMatrix[i-1];
          }

          for (j = 2; j < numOfPro ; ++j)
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
        else if (recvIndex == 15)
        {
          matrix[0][0]=recvupLeftNum;
        }
      }

      newMatrix[0][0]=((((matrix[0][0]*filter[0][0])+
                          (matrix[0][1]*filter[0][1])+
                          (matrix[0][2]*filter[0][2])+
                          (matrix[1][0]*filter[1][0])+
                          (matrix[1][1]*filter[1][1])+
                          (matrix[1][2]*filter[1][2])+
                          (matrix[2][0]*filter[2][0])+
                          (matrix[2][1]*filter[2][1])+
                          (matrix[2][2]*filter[2][2])))/normalize);

      newMatrix[0][numOfPro-1]=((((matrix[0][numOfPro-1]*filter[0][0])+
                                (matrix[0][numOfPro]*filter[0][1])+
                                (matrix[0][numOfPro+1]*filter[0][2])+
                                (matrix[1][numOfPro-1]*filter[1][0])+
                                (matrix[1][numOfPro]*filter[1][1])+
                                (matrix[1][numOfPro+1]*filter[1][2])+
                                (matrix[2][numOfPro-1]*filter[2][0])+
                                (matrix[2][numOfPro]*filter[2][1])+
                                (matrix[2][numOfPro+1]*filter[2][2])))/normalize);

      newMatrix[numOfPro-1][0]=((((matrix[numOfPro-1][0]*filter[0][0])+
                                  (matrix[numOfPro-1][1]*filter[0][1])+
                                  (matrix[numOfPro-1][2]*filter[0][2])+
                                  (matrix[numOfPro][0]*filter[1][0])+
                                  (matrix[numOfPro][1]*filter[1][1])+
                                  (matrix[numOfPro][2]*filter[1][2])+
                                  (matrix[numOfPro+1][0]*filter[2][0])+
                                  (matrix[numOfPro+1][1]*filter[2][1])+
                                  (matrix[numOfPro+1][2]*filter[2][2])))/normalize);

      newMatrix[numOfPro-1][numOfPro-1]=((((matrix[numOfPro-1][numOfPro-1]*filter[0][0])+
                                          (matrix[numOfPro-1][numOfPro]*filter[0][1])+
                                          (matrix[numOfPro-1][numOfPro+1]*filter[0][2])+
                                          (matrix[numOfPro][numOfPro-1]*filter[1][0])+
                                          (matrix[numOfPro][numOfPro]*filter[1][1])+
                                          (matrix[numOfPro][numOfPro+1]*filter[1][2])+
                                          (matrix[numOfPro+1][numOfPro-1]*filter[2][0])+
                                          (matrix[numOfPro+1][numOfPro]*filter[2][1])+
                                          (matrix[numOfPro+1][numOfPro+1]*filter[2][2])))/normalize);
      int sum=0;
      for (i = 0; i < numOfPro+2; i++)
      {
        for (j = 0; j < numOfPro+2; j++)
        {

          if (i==0 || i==numOfPro+1 || j==0 || j == numOfPro+1 )
          {
            matrix[i][j]= 255;
          }
          else{
            if(newMatrix[i-1][j-1]>=255)
            {
              matrix[i][j]=255;
            }else if (newMatrix[i-1][j-1]<=0)
            {
              matrix[i][j]=0;
            }
            else
            {
              matrix[i][j]=newMatrix[i-1][j-1];
            }  
            sum+=matrix[i][j];
          }
        }
      }
      int flag=sum/(numOfPro*numOfPro);
      if (((flag>=240 && flag<=255) || (flag<=15 && flag>=0) && foundS==0))
      {
        foundS = s;
      }
      /*int buff=-10;
      MPI_Isend(&buff, 1, MPI_INT, inbuf[UPLEFT], tag, MPI_COMM_WORLD, &newReqs[0]);
      MPI_Isend(&buff, 1, MPI_INT, inbuf[UP], tag, MPI_COMM_WORLD, &newReqs[1]);
      MPI_Isend(&buff, 1, MPI_INT, inbuf[UPRIGHT], tag, MPI_COMM_WORLD, &newReqs[2]);
      MPI_Isend(&buff, 1, MPI_INT, inbuf[LEFT], tag, MPI_COMM_WORLD, &newReqs[3]);
      MPI_Isend(&buff, 1, MPI_INT, inbuf[RIGHT], tag, MPI_COMM_WORLD, &newReqs[4]);
      MPI_Isend(&buff, 1, MPI_INT, inbuf[DOWNLEFT], tag, MPI_COMM_WORLD, &newReqs[5]);
      MPI_Isend(&buff, 1, MPI_INT, inbuf[DOWN], tag, MPI_COMM_WORLD, &newReqs[6]);
      MPI_Isend(&buff, 1, MPI_INT, inbuf[DOWNRIGHT], tag, MPI_COMM_WORLD, &newReqs[7]);
      
      MPI_Irecv(&buff, 1, MPI_INT, inbuf[DOWNRIGHT], tag, MPI_COMM_WORLD, &newReqs[8]);   
      MPI_Irecv(&buff, 1, MPI_INT, inbuf[DOWN], tag, MPI_COMM_WORLD, &newReqs[9]);
      MPI_Irecv(&buff, 1, MPI_INT, inbuf[DOWNLEFT], tag, MPI_COMM_WORLD, &newReqs[10]);
      MPI_Irecv(&buff, 1, MPI_INT, inbuf[RIGHT], tag, MPI_COMM_WORLD, &newReqs[11]);
      MPI_Irecv(&buff, 1, MPI_INT, inbuf[LEFT], tag, MPI_COMM_WORLD, &newReqs[12]);
      MPI_Irecv(&buff, 1, MPI_INT, inbuf[UPRIGHT], tag, MPI_COMM_WORLD, &newReqs[13]);
      MPI_Irecv(&buff, 1, MPI_INT, inbuf[UP], tag, MPI_COMM_WORLD, &newReqs[14]);
      MPI_Irecv(&buff, 1, MPI_INT, inbuf[UPLEFT], tag, MPI_COMM_WORLD, &newReqs[15]);

      MPI_Waitall(16, newReqs,newStats);*/
    }
//------------------------------------------------------------------------------------------------/

    //printf("###bbbbbbbbbbbb######\n" );
    //sleep(2);
    printf("Rank :::: %d s= %d     foundS=%d \n",rank ,s ,foundS);
    for (i = 0; i < numOfPro+2; ++i)
    {
      for (j = 0; j < numOfPro+2; ++j)
      {
        printf("%d  ", matrix[i][j]);
      }
      printf("\n");
    }
    printf("######################\n");
    for (i = 0; i < numOfPro; ++i)
    {
      for (j = 0; j < numOfPro; ++j)
      {
        printf("%d  ", newMatrix[i][j]);
      }
      printf("\n");
    }
    printf("miamoyRank :::: %d\n",rank);
  }
  else{
      printf("Must specify %d processors. Terminating.\n",SIZE);
   }
   
  MPI_Finalize();
}
