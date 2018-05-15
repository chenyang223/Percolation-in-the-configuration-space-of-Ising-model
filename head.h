#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <mpi.h>

using namespace std;

extern const int kNumSample;
extern const int kSingle;
extern const int kNumDistance;

extern int g_num_nearest_neighbor;
extern int g_ising_mode;
extern int g_length_ising;
extern int g_num_ising_site;

extern clock_t g_start_time,g_current_time;
extern int g_myrank;

extern int *g_ising_sites;     // Ising grids
extern int **g_sample;
extern int *g_sample_ptr;

extern double g_T;//Temperature
extern double g_J;//coefficient of Ising modle

//common.cpp
double RandDouble(int i);
int RandInt(int i);

//ising.cpp
void SampleConfiguration();
void MonteCarloWolff();

//Net_percolation.cpp
int CalculateDelta();
double CalculateDeltaMC();
