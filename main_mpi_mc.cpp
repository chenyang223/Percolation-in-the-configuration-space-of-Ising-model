#include "head.h"

double g_T;
double g_J = 1.0;

const int kNumSample = 10000;
const int kSingle = -1;
const int kNumDistance = (kNumSample * (kNumSample - 1)) / 2;

int g_length_ising;
//when lenthg_ising changes, remember to update g_num_ising_site
int g_num_ising_site;

int g_ising_mode; 
//# of the nearest neighbor for one ising site
//for normal ising model, the value is 4
//the value here can be 3,4 and 6
//the minus value of g_ising_mode means 3d ising model and for 3d case g_ising_mode is -6
//but 3d case haven't been tested in this program
int g_num_nearest_neighbor;//the abs value of g_ising_mode 

int *g_ising_sites;		//
int **g_sample;			// array to store samples
int *g_sample_ptr;
//Newman's method in fast percolation algorithms
//pointer showing how this node is connected to clusters
//positive: the value is the root or subroot(whose g_sample_ptr points to root)
//negative: opposite value of cluster size

clock_t g_start_time, g_current_time;
int npes, g_myrank;
//npes: number of cores, g_myrank: the rank of current process

/*
//available for Monte Carlo simulation to test whether the program is correct
double Tmin = 2.0;
double Tcurrent;
double delta_T = 0.05;
double Tmax;
double m;
double H;
*/

//allocate memory g_ising_sites, g_sample and g_sample_ptr
void AllocateMemoryToArray()
{
	g_ising_sites = new int[g_num_ising_site];
	g_sample = new int *[kNumSample];
	for (int i = 0; i < kNumSample; i++)
		g_sample[i] = new int[g_num_ising_site];
	g_sample_ptr = new int[kNumSample];
}

void DeleteMemoryToArray()
{
	delete[] g_ising_sites;
	for (int i = 0; i < kNumSample; i++)
		delete[] g_sample[i];
	delete[] g_sample;
	delete[] g_sample_ptr;
}

int main()
{

	//MPI
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &g_myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);

	//change parameters at this part
	int repeat_time = 1;
	g_ising_mode = 6;
	g_num_nearest_neighbor = abs(g_ising_mode);

	int length_T_list = 20;
	float *T_list = new float[length_T_list];
	for (int i = 0; i < length_T_list; i++)
		T_list[i] = 3.55 + 0.01 * i;
	int length_of_length_list = 3;
	int *length_list = new int[length_of_length_list];
	for (int i = 0; i < length_of_length_list; i++)
		length_list[i] = pow(2, 3 + i);


	double difference_current = 0.0; // save the value of \Delta temporarily
	int length_diff = 4;
	double *difference_biggest_clu = new double[length_diff]; //array store the 1st to 4th power of diff.
	double **difference_list = new double *[length_diff];
	for (int i = 0; i < length_diff; i++)
	{
		difference_biggest_clu[i] = 0.0;
		difference_list[i] = new double[length_T_list];
	}

	srand((unsigned)time(NULL) + g_myrank);

	for (int i = 0; i < length_of_length_list; i++)
	{
		g_length_ising = length_list[i];
		if(g_ising_mode > 0)
			g_num_ising_site = g_length_ising * g_length_ising;
		else
			g_num_ising_site = g_length_ising * g_length_ising * g_length_ising;

		for (int j = 0; j < length_T_list; j++)
		{
			g_start_time = clock();
			g_T = T_list[j];
			for (int k = 0; k < length_diff; k++)
				difference_biggest_clu[k] = 0.0;

			for (int k = 0; k < repeat_time; k++)
			{
				if (g_myrank == 0)
				{
					cout << g_length_ising << "    " << g_T;
					cout << "     " << k << "/" << repeat_time << endl;
				}
				AllocateMemoryToArray();
				SampleConfiguration();
				if (g_myrank == 0)
				{
					g_current_time = clock();
					cout << "g_sample finished-----process" << g_myrank << "----time:";
					cout << (double)(g_current_time - g_start_time) / CLOCKS_PER_SEC << endl;
				}
				difference_current = CalculateDeltaMC();
				difference_current /= kNumSample;
				for (int l = 0; l < length_diff; l++)
					difference_biggest_clu[l] += pow(difference_current, l + 1);
				DeleteMemoryToArray();
			}
			for (int k = 0; k < length_diff; k++)
				difference_biggest_clu[k] /= (repeat_time);

			MPI_Barrier(MPI_COMM_WORLD);
			double **difference_of_biggest_clu_set = new double *[length_diff];
			if (g_myrank == 0)
			{
				for (int k = 0; k < length_diff; k++)
					difference_of_biggest_clu_set[k] = new double[npes];
			}
			for (int k = 0; k < length_diff; k++)
				MPI_Gather(&difference_biggest_clu[k], 1, MPI_DOUBLE,
						   difference_of_biggest_clu_set[k], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if (g_myrank == 0)
			{
				for (int k = 0; k < length_diff; k++)
				{
					difference_biggest_clu[k] = 0;
					for (int l = 0; l < npes; l++)
					{
						difference_biggest_clu[k] += difference_of_biggest_clu_set[k][l];
					}
					difference_biggest_clu[k] /= npes;
					difference_list[k][j] = difference_biggest_clu[k];
					delete [] difference_of_biggest_clu_set[k];
				}
			}
			delete [] difference_of_biggest_clu_set;
		}
		//file output
		if (g_myrank == 0)
		{
			stringstream file_name_stream;
			string filename_difference_biggest_clu;
			file_name_stream << g_length_ising << "_" << g_ising_mode << "_";
			file_name_stream << "difference_of_biggest_cluster_MC.txt";
			file_name_stream >> filename_difference_biggest_clu;
			ofstream outfile_difference_biggest_clu;
			outfile_difference_biggest_clu.open(filename_difference_biggest_clu.c_str());
			for (int j = 0; j < length_T_list; j++)
			{
				outfile_difference_biggest_clu << T_list[j] << " ";
				for (int k = 0; k < length_diff; k++)
					outfile_difference_biggest_clu << difference_list[k][j] << " ";
				outfile_difference_biggest_clu << endl;
			}
		}
	}
	for (int i = 0; i < length_diff; i++)
	{
		delete[] difference_list[i];
	}
	delete[] difference_list;
	delete[] difference_biggest_clu;
	delete[] T_list;
	delete[] length_list;
	MPI_Finalize();
	return 0;
}