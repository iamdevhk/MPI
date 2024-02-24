#include <iostream>
#include "Timer.h"
#include "mpi.h"
#include <stdlib.h> // atoi
#include <stdio.h>
#include <math.h>
#include <omp.h>

int default_size = 100; // the default system size
int defaultCellWidth = 8;
double c = 1.0;  // wave speed
double dt = 0.1; // time quantum
double dd = 2.0; // change in system

using namespace std;

int main(int argc, char *argv[])
{
  int mpi_rank = 0; // used by MPI
  int mpi_size;     // number of processors

  // verify arguments
  if (argc != 5)
  {
    cerr << "usage: Wave2D size max_time interval" << endl;
    return -1;
  }

  int size = atoi(argv[1]);
  int max_time = atoi(argv[2]);
  int interval = atoi(argv[3]);
  int threadCount = atoi(argv[4]);

  if (size < 100 || max_time < 3 || interval < 0)
  {
    cerr << "usage: Wave2D size max_time interval" << endl;
    cerr << "       where size >= 100 && time >= 3 && interval >= 0" << endl;
    return -1;
  }

  //initialize openMP and MPI
  omp_set_num_threads(threadCount);
  MPI_Init(&argc, &argv); // start MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Status status;

  // create a simulation space
  double z[3][size][size];
  for (int p = 0; p < 3; p++)
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++)
        z[p][i][j] = 0.0; // no wave

  // start a timer
  Timer time;
  time.start();

  // time = 0;
  // initialize the simulation space --- z[0][][]
  int weight = size / default_size;
  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      if (i > 40 * weight && i < 60 * weight && j > 40 * weight && j < 60 * weight)
      {
        z[0][i][j] = 20.0;
      }
      else
      {
        z[0][i][j] = 0.0;
      }
    }
  }

// time = 1
// calculate z[1][][] --- cells not on edge
// IMPLEMENT BY YOURSELF !!!
#pragma omp parallel for
  for (int xItr = 1; xItr < size - 1; xItr++) //iterate row x
  {
    for (int yItr = 1; yItr < size - 1; yItr++) //iterate column y
    {
      z[1][xItr][yItr] = z[0][xItr][yItr] + (c * c) / 2 * ((dt * dt) / (dd * dd)) * (z[0][xItr + 1][yItr] + z[0][xItr - 1][yItr] + z[0][xItr][yItr + 1] + z[0][xItr][yItr - 1] - 4.0 * z[0][xItr][yItr]);
    }
  }

  int partitionedStrip = size / mpi_size; // partitioned stripe
  int stripesArr[mpi_size];
  int remainder = size % mpi_size; //leftover which cannot be evenly distributed to mpi process

  // add values to stripesArr based on remainder
  for (int itr = 0; itr < mpi_size; itr++)
  {
    if (itr < remainder)
    {
      stripesArr[itr] = partitionedStrip + 1;
    }
    else
    {
      stripesArr[itr] = partitionedStrip;
    }
  }

  // simulate wave diffusion from time = 2
  // IMPLEMENT BY YOURSELF !!!
  for (int timeIterator = 2; timeIterator < max_time; timeIterator++)
  {
    int curWave = timeIterator % 3;        // current wave
    int prevWave = (timeIterator + 1) % 3; // previous wave
    int priWave = (timeIterator + 2) % 3;  // previous to previous - prior wave

    // for the remainder case:
    if (mpi_rank == 0)
    { // master exchanges only the right data to rank+1                  
      if (mpi_size > 1) // checks if more than 1 processors exist
      {

        // sends data of last lower boundary
        MPI_Send(*(*(z + priWave) + stripesArr[mpi_rank] - 1), size, MPI_DOUBLE, mpi_rank + 1, 0, MPI_COMM_WORLD);
        // recievies data for lower boundary+1
        MPI_Recv(*(*(z + priWave) + stripesArr[mpi_rank]), size, MPI_DOUBLE, mpi_rank + 1, 0, MPI_COMM_WORLD, &status);
      }
    }
    // last processor exchanges its left boundary data to rank-1
    else if (mpi_rank == mpi_size - 1)
    {                   
      if (mpi_size > 1) 
      {

        // sends the data of upper first boundary
        MPI_Send(*(*(z + priWave) + (mpi_rank * stripesArr[mpi_rank]) + remainder), size, MPI_DOUBLE, mpi_rank - 1, 0, MPI_COMM_WORLD);
        // recieves the data for upper boundary-1
        MPI_Recv(*(*(z + priWave) + (mpi_rank * stripesArr[mpi_rank]) + remainder - 1), size, MPI_DOUBLE, mpi_rank - 1, 0, MPI_COMM_WORLD, &status);
      }
    }
    else
    {                   // other processors exchange their left and right boundary data to rank-1 and rank+1
      if (mpi_size > 1) 
      {
        if (mpi_rank < remainder)
        {
          // sends data of first upper boundary 
          MPI_Send(*(*(z + priWave) + stripesArr[mpi_rank] * mpi_rank), size, MPI_DOUBLE, mpi_rank - 1, 0, MPI_COMM_WORLD);
          // recieves data for upper boundary-1
          MPI_Recv(*(*(z + priWave) + (stripesArr[mpi_rank] * mpi_rank) - 1), size, MPI_DOUBLE, mpi_rank - 1, 0, MPI_COMM_WORLD, &status);
          // sends data of last lower boundary 
          MPI_Send(*(*(z + priWave) + (stripesArr[mpi_rank] * (mpi_rank + 1)) - 1), size, MPI_DOUBLE, mpi_rank + 1, 0, MPI_COMM_WORLD);
          // recieves data for lower boundary+1
          MPI_Recv(*(*(z + priWave) + stripesArr[mpi_rank] * (mpi_rank + 1)), size, MPI_DOUBLE, mpi_rank + 1, 0, MPI_COMM_WORLD, &status);
        }
        else
        {
          // sends data of first upper boundary
          MPI_Send(*(*(z + priWave) + (stripesArr[mpi_rank] * mpi_rank) + remainder), size, MPI_DOUBLE, mpi_rank - 1, 0, MPI_COMM_WORLD);
          // recieves data  upper boundary-1
          MPI_Recv(*(*(z + priWave) + (stripesArr[mpi_rank] * mpi_rank) + remainder - 1), size, MPI_DOUBLE, mpi_rank - 1, 0, MPI_COMM_WORLD, &status);
          // sends data of last lower boundary 
          MPI_Send(*(*(z + priWave) + (stripesArr[mpi_rank] * (mpi_rank + 1)) + remainder - 1), size, MPI_DOUBLE, mpi_rank + 1, 0, MPI_COMM_WORLD);
          // recieves data for lower boundary+1
          MPI_Recv(*(*(z + priWave) + (stripesArr[mpi_rank] * (mpi_rank + 1)) + remainder), size, MPI_DOUBLE, mpi_rank + 1, 0, MPI_COMM_WORLD, &status);
        }
      }
    }

    // computes z data for current timeIterator-- using parallelization to have multiple threads wokring on different parts of the grid

    if (mpi_rank < remainder)
    {
#pragma omp parallel for
      for (int xItr = (mpi_rank * stripesArr[mpi_rank]); xItr < (mpi_rank + 1) * stripesArr[mpi_rank]; xItr++) //iterate rows x
      {
        if (xItr == 0 || xItr == size - 1) // skip boundary rows
        {
          continue;
        }
        for (int yItr = 1; yItr < size - 1; yItr++) //iterate columns y
        {
          z[curWave][xItr][yItr] = 2.0 * z[priWave][xItr][yItr] - z[prevWave][xItr][yItr] + c * c * dt * dt / (dd * dd) * (z[priWave][xItr + 1][yItr] + z[priWave][xItr - 1][yItr] + z[priWave][xItr][yItr + 1] + z[priWave][xItr][yItr - 1] - 4.0 * z[priWave][xItr][yItr]);
        }
      }
    }
    else // no additonal row
    {
#pragma omp parallel for
      for (int xItr = ((mpi_rank * stripesArr[mpi_rank]) + remainder); xItr < (((mpi_rank + 1) * stripesArr[mpi_rank]) + remainder); xItr++) //iterate row x
      {
        if (xItr == 0 || xItr == size - 1) // skip boundary rows
        {
          continue;
        }
        for (int yItr = 1; yItr < size - 1; yItr++) //iterate column y
        {
          z[curWave][xItr][yItr] = 2.0 * z[priWave][xItr][yItr] - z[prevWave][xItr][yItr] + c * c * dt * dt / (dd * dd) * (z[priWave][xItr + 1][yItr] + z[priWave][xItr - 1][yItr] + z[priWave][xItr][yItr + 1] + z[priWave][xItr][yItr - 1] - 4.0 * z[priWave][xItr][yItr]);
        }
      }
    }

    // printing the data for the interval
    if (interval != 0 && timeIterator % interval == 0)
    {

      // send data back to master from the workers
      if (mpi_rank != 0)
      {
        if (mpi_rank < remainder)
        {
          MPI_Send(*(*(z + curWave) + (stripesArr[mpi_rank] * mpi_rank)), stripesArr[mpi_rank] * size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        else
        {
          MPI_Send(*(*(z + curWave) + (stripesArr[mpi_rank] * mpi_rank) + remainder), stripesArr[mpi_rank] * size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
      }
      else
      {
        // master recieves all the data from workers
        for (int rank = 1; rank < mpi_size; rank++)
        {
          if (rank < remainder)
          {
            MPI_Recv(*(*(z + curWave) + rank * stripesArr[rank]), stripesArr[rank] * size, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, &status);
          }
          else
          {
            MPI_Recv(*(*(z + curWave) + rank * stripesArr[rank] + remainder), stripesArr[rank] * size, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, &status);
          }
        }
        // printing z data for current timeIterator
        cout << timeIterator << endl;
        for (int xItr = 0; xItr < size; xItr++) //iterate rows x
        {
          for (int yItr = 0; yItr < size; yItr++) //iterate columns y
          {
            cout << z[curWave][yItr][xItr] << " ";
          }
          cout << endl;
        }
        cout << endl;
      }
    }
  }
  // end of simulation

  MPI_Finalize(); // closes all mpi related resources used 

  
  if (mpi_rank == 0)
  {
    //print for each rank
    for (int rankItr = 0; rankItr < mpi_size; rankItr++)
    {
      if (rankItr < remainder)
      {
        cerr << "rank[" << rankItr << "]'s range " << rankItr * stripesArr[rankItr] << "~" << (rankItr + 1) * stripesArr[rankItr] - 1 << endl;
      }
      else
      {
        cerr << "rank[" << rankItr << "]'s range " << rankItr * stripesArr[rankItr] + remainder << "~" << (rankItr + 1) * stripesArr[rankItr] + (remainder - 1) << endl;
      }
    }
    cerr << "Elapsed time = " << time.lap() << endl; // prints the elapsed time
  }
  return 0;
}
