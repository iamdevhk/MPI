#include <iostream>
#include "Timer.h"
#include <stdlib.h> // atoi
#include <stdio.h>
#include <math.h>

int default_size = 100; // the default system size
int defaultCellWidth = 8;
double c = 1.0;  // wave speed
double dt = 0.1; // time quantum
double dd = 2.0; // change in system

using namespace std;

int main(int argc, char *argv[])
{
  // verify arguments
  if (argc != 4)
  {
    cerr << "usage: Wave2D size max_time interval" << endl;
    return -1;
  }

  int size = atoi(argv[1]);
  int max_time = atoi(argv[2]);
  int interval = atoi(argv[3]);

  if (size < 100 || max_time < 3 || interval < 0)
  {
    cerr << "usage: Wave2D size max_time interval" << endl;
    cerr << "       where size >= 100 && time >= 3 && interval >= 0" << endl;
    return -1;
  }

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
  // initialize the simulation space: calculate z[0][][]
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
  // calculate z[1][][]
  // cells not on edge
  // IMPLEMENT BY YOURSELF !!!
  for (int xItr = 1; xItr < size - 1; xItr++) // iterate rows x
  {
    for (int yItr = 1; yItr < size - 1; yItr++) // iterate columns y
    {
      // updates the wave value at xItr,yItr using the wave equation
      z[1][xItr][yItr] = z[0][xItr][yItr] + (c * c) / 2 * ((dt * dt) / (dd * dd)) * (z[0][xItr + 1][yItr] + z[0][xItr - 1][yItr] + z[0][xItr][yItr + 1] + z[0][xItr][yItr - 1] - 4.0 * z[0][xItr][yItr]);
    }
  }

  // simulate wave diffusion from time = 2
  // IMPLEMENT BY YOURSELF !!!
  // iterates from 2 to max_time entered by the user
  for (int timeIterator = 2; timeIterator < max_time; timeIterator++)
  {

    int curWave = timeIterator % 3;        // current wave
    int prevWave = (timeIterator + 1) % 3; // previous wave
    int priWave = (timeIterator + 2) % 3;  // previous to previous - prior wave
    // iterate over the grid cells
    for (int xItr = 1; xItr < size - 1; xItr++) // iterate rows x
    {
      for (int yItr = 1; yItr < size - 1; yItr++) // iterate columns y
      {
        // updates the wave value of xItr - yItr using the given wave equation
        z[curWave][xItr][yItr] = 2.0 * z[priWave][xItr][yItr] - z[prevWave][xItr][yItr] + c * c * dt * dt / (dd * dd) * (z[priWave][xItr + 1][yItr] + z[priWave][xItr - 1][yItr] + z[priWave][xItr][yItr + 1] + z[priWave][xItr][yItr - 1] - 4.0 * z[priWave][xItr][yItr]);
      }
    }
    if (interval != 0 && timeIterator % interval == 0)
    {
      cout << timeIterator << endl;
      for (int xItr = 0; xItr < size; xItr++) // iterate rows x
      {
        for (int yItr = 0; yItr < size; yItr++) // iterate columns y
        {
          cout << z[curWave][yItr][xItr] << " ";
        }
        cout << endl;
      }
      cout << endl;
    }
  }
  // end of simulation

  // finish the timer
  cerr << "Elapsed time = " << time.lap() << endl;
  return 0;
}
