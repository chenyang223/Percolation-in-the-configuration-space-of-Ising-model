//all the abbreviation "clu" in this file means cluster
#include "head.h"

typedef struct DistanceMatrixElement
{
	int index1;
	int index2; // the labels of the vector being multiplied
	int dist;   // the opposite number of innerproduct between g_sample i and j
} DistanceMatrixElement;

/* //available for calculating center of cluster and effective radius
//If wanting to activate this, remember to activate the memory allocation
// and deletion and computation part in all function
double **center_of_clu;
double *effective_radius;
*/

int mc_bond_coefficient = 20;
int num_mc_bonds = mc_bond_coefficient * kNumSample;
//available for monte carlo network construction
//the # of bonds randomly generated is mc_bond_coefficients*g_num_ising_site
int num_add_bonds = 0;
//the # of bonds that is going to be added

int *num_of_sizes;
//array to save the number of clusters of different sizes
//the index means the size of the cluster
//num_of_sizes[x] means the # of cluster of size x
DistanceMatrixElement *distances_sample;

int big_clu_size[5] = {1, 1, 1, 1, 1}; //sizes of the 1st to 5th biggest cluster
int big_clu_root[2] = {0, 1};		   //root for the 1st and 2nd biggest cluster

static void AllocateMemoryToArray()
{
	num_of_sizes = new int[kNumSample + 1];
	distances_sample = new DistanceMatrixElement[kNumDistance];
	/*
	center_of_clu = new double *[kNumSample];
	for(int i = 0; i < kNumSample; i++)
		center_of_clu[i] = new double [g_num_ising_site];
	effective_radius = new double [kNumSample];
	*/
}

static void AllocateMemoryToArrayMC()
{
	num_of_sizes = new int[kNumSample + 1];
	distances_sample = new DistanceMatrixElement[num_mc_bonds];
}

static void DeleteMemoryToArray()
{
	delete[] num_of_sizes;
	delete[] distances_sample;
	/*
	for(int i = 0; i < kNumSample; i++)
		delete [] center_of_clu[i];
	delete [] center_of_clu;
	delete [] effective_radius;
	*/
}

//not only findroot, but also path compression!!
int FindRoot(int i)
{
	if (g_sample_ptr[i] < 0)
		return i;
	return g_sample_ptr[i] = FindRoot(g_sample_ptr[i]);
}

void PrintDistance()
{
	string file_name_distance = "distance.txt";
	stringstream file_name_stream;
	file_name_stream << setprecision(3) << g_T << "_";
	file_name_stream << g_length_ising << "_" << file_name_distance;
	file_name_stream >> file_name_distance;
	ofstream out_distance;
	out_distance.open(file_name_distance.c_str());

	string file_name_distace_average = "distance_avergae.txt";
	file_name_stream.str("");
	file_name_stream << setprecision(3) << g_T << "_";
	file_name_stream << g_length_ising << "_" << file_name_distace_average;
	file_name_stream >> file_name_distace_average;
	ofstream out_distance_avg;
	out_distance.open(file_name_distace_average.c_str());

	double distance_sum = 0;

	for (int i = 0; i < kNumDistance; i++)
	{
		out_distance << distances_sample[i].dist << " ";
		distance_sum += distances_sample[i].dist;
		out_distance_avg << distance_sum / (i + 1) << " ";
		//printf("%d %I64d\n", distances_sample[i].prod, distance_sum / (i + 1));

		out_distance.close();
		out_distance_avg.close();
	}
}

//initialize g_sample_ptr, num_of_sizes, distances_sample and num_add_bonds
//There are three Initialize function in this file
//1. Initialize: the original initialzing function which calculating the distance directly
//2. InitializeOptimized: it calculate the distance by xor
//3. InitializeMC: it calculate the distance by xor and is available for MC case
void Initialize()
{
	/*
	FILE *fp1;
	FILE *fp2;
	char name11[] = "+distance.txt";
	char name22[] = "-distance.txt";
	char name1[30];
	char name2[30];
	sprintf(name1, "%.3lf%s", g_T, name11);
	sprintf(name2, "%.3lf%s", g_T, name22);
	fp1 = fopen(name1, "w+");
	fp2 = fopen(name2, "w+");
	*/

	num_add_bonds = kNumSample * (kNumSample - 1) / 2;
	DistanceMatrixElement *pp_temp = distances_sample;
	g_sample_ptr[0] = kSingle;
	//effective_radius[0] = 0.0;
	num_of_sizes[0] = 0;

	//for (int j = 0; j < g_num_ising_site; j++)
	//	center_of_clu[0][j] = g_sample[0][j];

	for (int i = 1; i < kNumSample; i++)
	{
		g_sample_ptr[i] = kSingle;
		//effective_radius[i] = 0.0;
		num_of_sizes[i] = 0;

		for (int j = 0; j < kNumSample; j++)
		{
			//if (j < g_num_ising_site)
			//	center_of_clu[i][j] = g_sample[i][j];
			if (j < i)
			{
				pp_temp->index1 = i;
				pp_temp->index2 = j;
				pp_temp->dist = 0;
				for (int k = 0; k < g_num_ising_site; k++)
					pp_temp->dist -= g_sample[i][k] * g_sample[j][k];

				/*
				if (pp_temp->dist < 0)
					fprintf(fp2, "%d ", pp_temp->dist);
				else
					fprintf(fp1, "%d ", pp_temp->dist);
				*/
				pp_temp++;
			}
		}
	}
	num_of_sizes[kNumSample] = 0;
	num_of_sizes[1] = kNumSample;
}

//use xor to calculate the innerproduct of samples
//the size of the biggest data strcture unsigned long long int is 64 bit
//seperate the g_sample into many fragments of size 64
//calculate xor and count the # of '1' in the result
void InitializeOptimized()
{

	//cout << "---------------Original Initialization------------" << endl;
	//clock_t time1 = clock();
	DistanceMatrixElement *pp_temp = distances_sample;
	num_add_bonds = kNumSample * (kNumSample - 1) / 2;

	int size_of_fragment = 64; //size of unsigned long long int
	int num_sample_fragment = g_num_ising_site / size_of_fragment;
	int remainder_sample_fragment = g_num_ising_site % size_of_fragment;
	if (num_sample_fragment == 0 || remainder_sample_fragment != 0)
		num_sample_fragment++;
	unsigned long long **sample_binary;
	//sample_binary is an array saving samples as binary number
	//1 and 0 means two kinds of spin
	sample_binary = new unsigned long long *[kNumSample];
	//transform from g_sample to sample_binary
	for (int i = 0; i < kNumSample; i++)
	{
		g_sample_ptr[i] = kSingle;
		num_of_sizes[i] = 0;
		sample_binary[i] = new unsigned long long[num_sample_fragment];
		sample_binary[i][num_sample_fragment - 1] = 0;
		//this step is to guarantee the bit not containing ising site is 0
		for (int j = 0; j < g_num_ising_site; j++)
		{
			int index = j / size_of_fragment;
			sample_binary[i][index] = sample_binary[i][index] << 1;
			if (g_sample[i][j] == 1)
				sample_binary[i][index]++;
		}
	}
	num_of_sizes[kNumSample] = 0;
	num_of_sizes[1] = kNumSample;

	//clock_t time2 = clock();

	for (int i = 1; i < kNumSample; i++)
	{
		for (int j = 0; j < i; j++)
		{
			pp_temp->index1 = i;
			pp_temp->index2 = j;
			pp_temp->dist = -remainder_sample_fragment;
			for (int k = 0; k < num_sample_fragment; k++)
			{
				unsigned long long xor_result =
					sample_binary[i][k] ^ sample_binary[j][k];
				int innerprod = 0;
				// counting the number of "1" in xor_result
				while (xor_result != 0)
				{
					xor_result = xor_result & (xor_result - 1);
					innerprod++;
				}
				pp_temp->dist += innerprod;
			}
			pp_temp++;
		}
	}
	for (int i = 0; i < kNumSample; i++)
	{
		delete[] sample_binary[i];
	}
	delete[] sample_binary;
	//clock_t time3 = clock();
	//cout << endl;
	//cout << double(time2 - time1) / CLOCKS_PER_SEC << " ";
	//cout << double(time3 - time2) / CLOCKS_PER_SEC << endl;
}


void InitializeMC()
{
	//cout << "--------------------MC Initialization-------------" << endl;
	//clock_t time1 = clock();
	DistanceMatrixElement *pp_temp = distances_sample;
	num_add_bonds = num_mc_bonds;

	int size_of_fragment = 64; //size of unsigned long long int
	int num_sample_fragment = g_num_ising_site / size_of_fragment;
	int remainder_sample_fragment = g_num_ising_site % size_of_fragment;
	if (num_sample_fragment == 0 || remainder_sample_fragment != 0)
		num_sample_fragment++;
	unsigned long long **sample_binary;
	sample_binary = new unsigned long long *[kNumSample];

	for (int i = 0; i < kNumSample; i++)
	{
		g_sample_ptr[i] = kSingle;
		num_of_sizes[i] = 0;
		sample_binary[i] = new unsigned long long[num_sample_fragment];
		sample_binary[i][num_sample_fragment - 1] = 0;
		//this step is to guarantee the bit not containing ising site is 0
		for (int j = 0; j < g_num_ising_site; j++)
		{
			int index = j / size_of_fragment;
			sample_binary[i][index] = sample_binary[i][index] << 1;
			if (g_sample[i][j] == 1)
				sample_binary[i][index]++;
		}
	}
	num_of_sizes[kNumSample] = 0;
	num_of_sizes[1] = kNumSample;

	//clock_t time2 = clock();
	int addbond_index1;
	int addbond_index2;
	for (int i = 0; i < num_mc_bonds; i++)
	{
		addbond_index1 = RandInt(kNumSample);
		addbond_index2 = RandInt(kNumSample);
		if (addbond_index1 == addbond_index2)
		{
			i--;
			continue;
		}
		pp_temp->index1 = addbond_index1;
		pp_temp->index2 = addbond_index2;
		pp_temp->dist = -remainder_sample_fragment;
		for (int k = 0; k < num_sample_fragment; k++)
		{
			unsigned long long xor_result =
				sample_binary[addbond_index1][k] ^ sample_binary[addbond_index2][k];
			int innerprod = 0;
			// counting the number of "1" in xor_result
			while (xor_result != 0)
			{
				xor_result = xor_result & (xor_result - 1);
				innerprod++;
			}
			pp_temp->dist += innerprod;
		}
		pp_temp++;
	}
	for (int i = 0; i < kNumSample; i++)
	{
		delete[] sample_binary[i];
	}
	delete[] sample_binary;
	//clock_t time3 = clock();
	//cout << endl;
	//cout << "step 1 time:" << double(time2 - time1) / CLOCKS_PER_SEC <<" ";
	//cout << "step 2 time:" << double(time3 - time2) / CLOCKS_PER_SEC << endl;
}

//select the mid value as the pivot for Quicksort
//meanwhile make sure the order is that right>=left>=mid
DistanceMatrixElement SelectPivot(DistanceMatrixElement *a, int left, int right)
{
	int mid = left + (right - left) / 2;
	if (a[mid].dist > a[right].dist)
		swap(a[mid], a[right]);
	if (a[left].dist > a[right].dist)
		swap(a[left], a[right]);
	if (a[mid].dist > a[left].dist)
		swap(a[mid], a[left]);
	return a[left];
}

//special optimization to repeating elements:
//find out all repeating elements and put them on the two sides.
//then Swap them to tha adjacent places of the pivot
void QuicksortDist(DistanceMatrixElement *arr, int left, int right)
{
	typedef struct stack
	{
		int left;
		int right;
	} stack;
	//artifitial stack

	stack *stk_original = new stack[(2 * g_num_ising_site)];
	if (!stk_original)
	{
		printf("The memory is full\n");
		exit(1);
	}
	stack *stk = stk_original;
	int stk_ptr = 0;
	stk[stk_ptr].left = left;
	stk[stk_ptr].right = right;

	for (;;)
	{ //partition part
		if (stk_ptr < 0)
			break;

		int a, b; //pointer to the two sides
		a = stk[stk_ptr].left;
		b = stk[stk_ptr--].right;
		if (a >= b)
			continue;
		int i = a, j = b; //pointer for QuicksortDist
		int m = a, n = b; //pointer for the repeat elements
		int len_right = 0, len_left = 0;

		DistanceMatrixElement pivot = SelectPivot(arr, a, b);

		while (i != j)
		{
			while (arr[j].dist >= pivot.dist && j > i)
			{
				if (arr[j].dist == pivot.dist)
				{
					swap(arr[j], arr[n]);
					n--;
					len_right++;
				}
				j--;
			}
			if (j > i)
				arr[i++] = arr[j];
			while (arr[i].dist <= pivot.dist && j > i)
			{
				if (arr[i].dist == pivot.dist)
				{
					swap(arr[i], arr[m]);
					m++;
					len_left++;
				}
				i++;
			}
			if (j > i)
				arr[j--] = arr[i];
		}
		arr[i] = pivot;

		m = a;	 //leftest side
		n = i - 1; //adjacent to pivot
		for (int k = 0; k < len_left; k++)
			swap(arr[m++], arr[n--]);

		m = i + 1; //adjacent to pivot
		n = b;	 //rightest side
		for (int k = 0; k < len_right; k++)
			swap(arr[m++], arr[n--]);

		//new range to be sorted
		stk[++stk_ptr].left = a;
		stk[stk_ptr].right = i - 1 - len_left;
		stk[++stk_ptr].left = i + 1 + len_right;
		stk[stk_ptr].right = b;
	}
	delete[] stk_original;
}

//require that cluster root2 is joined into cluster root1
//available befor changing g_sample_ptr
//For ising, R_eff = g_num_ising_site*\sum_{i=1}^{n}(1-\bar{x}_i^2)
//g_num_ising_site: # of sites in the cluster, n: # of elments for the vector
/*
void UpdateCenterOfMassAndEffectiveRadius(int root1, int root2) {
	int num_clu1 = -g_sample_ptr[root1];
	int num_clu2 = -g_sample_ptr[root2];

	effective_radius[root1] = 0;
	for (int i = 0; i < g_num_ising_site; i++)
		center_of_clu[root1][i] = (center_of_clu[root1][i] * num_clu1\
			+ center_of_clu[root2][i]) / (num_clu1 + num_clu2);
	for (int i = 0; i < g_num_ising_site; i++)
		effective_radius[root1] += 1 - center_of_clu[root1][i]*\
			center_of_clu[root1][i];
	effective_radius[root1] = sqrt(effective_radius[root1]);
}
*/

//available before changing g_sample_ptr
//root1,root2 : the root for the two combining cluster
//num_bond_added : x - 1, where x is the # of bonds that has already been
//added, because the loop parameter stars from 0 not 1.
void UpdateNumOfSizesAndBigCluster(int root1, int root2, int num_bond_added)
{
	num_of_sizes[-g_sample_ptr[root1]]--;
	num_of_sizes[-g_sample_ptr[root2]]--;
	num_of_sizes[-(g_sample_ptr[root1] + g_sample_ptr[root2])]++;

	//calculate array big_cluster_size
	//scan from the possible biggest cluster
	int max_possible = min(num_bond_added + 2, kNumSample);
	int n = 0;
	for (int j = 0; j < max_possible; j++)
	{
		int k = max_possible - j;
		if (num_of_sizes[k])
		{
			if (n + num_of_sizes[k] <= 5)
			{
				for (int l = 0; l < num_of_sizes[k]; l++)
					big_clu_size[l + n] = k;
			}
			else
			{
				for (int l = 0; l < (5 - n); l++)
					big_clu_size[l + n] = k;
			}
			n += num_of_sizes[k];
		}
		if (n >= 5)
			break;
	}
	if (n < 5)
	{
		for (; n < 5; n++)
			big_clu_size[n] = 0;
	}
}

//work after changing g_sample_ptr and updating big_clu_size
//only available for length-2 big_clu_root
void UpdateBigClusterRoot(int root_new)
{
	if (g_sample_ptr[big_clu_root[0]] < 0 && g_sample_ptr[big_clu_root[1]] < 0)
	{
		if (root_new == big_clu_root[0])
			return;
		if (root_new == big_clu_root[1] &&
			g_sample_ptr[big_clu_root[1]] < g_sample_ptr[big_clu_root[0]])
		{
			swap(big_clu_root[1], big_clu_root[0]);
			return;
		}
		if (g_sample_ptr[root_new] >= g_sample_ptr[big_clu_root[1]])
			return;
		if (g_sample_ptr[root_new] >= g_sample_ptr[big_clu_root[0]])
			big_clu_root[1] = root_new;
		else
		{
			big_clu_root[1] = big_clu_root[0];
			big_clu_root[0] = root_new;
		}
	}
	else
	{
		int root_retained;
		if (g_sample_ptr[big_clu_root[0]] < 0)
		{
			root_retained = big_clu_root[0];
		}
		else
		{
			root_retained = big_clu_root[1];
		}
		if (g_sample_ptr[root_retained] == -big_clu_size[0])
		{
			big_clu_root[0] = root_retained;
			for (int i = 0; i < kNumSample; i++)
			{
				if (g_sample_ptr[FindRoot(i)] == -big_clu_size[1])
				{
					big_clu_root[1] = FindRoot(i);
					break;
				}
			}
		}
		else
		{
			big_clu_root[1] = root_retained;
			for (int i = 0; i < kNumSample; i++)
			{
				if (g_sample_ptr[FindRoot(i)] == -big_clu_size[0])
				{
					big_clu_root[0] = FindRoot(i);
					break;
				}
			}
		}
	}
}

//method following Newman's paper about fast percolation algorithms
int BuildNetwork()
{
	int site1, site2;
	int root1, root2;
	int num_clu = kNumSample;
	double distance_between_two_biggest_clu = 0;
	int temp = big_clu_size[1];
	double temp1 = distance_between_two_biggest_clu;
	int difference_of_biggest_clu = 0; //
	int size_of_old_biggest_clu;

	/*
    string file_name_clu_dis;
    file_name_stream << setprecision(3) << g_T << "_" << g_length_ising << "_";
    file_name_stream << "#1,r1,#2,r2,dis.txt";
    file_name_stream >> file_name_clu_dis;
    ofstream outfile_clu_dis;
    outfile_clu_dis.open(file_name_clu_dis.c_str());
    */
	//adding bonds
	for (int i = 0; i < num_add_bonds; i++)
	{
		site1 = distances_sample[i].index1;
		site2 = distances_sample[i].index2;
		root1 = FindRoot(site1);
		root2 = FindRoot(site2);

		//combining clusters
		if (root1 != root2)
		{
			size_of_old_biggest_clu = big_clu_size[0];
			num_clu--;
			UpdateNumOfSizesAndBigCluster(root1, root2, i);
			if (g_sample_ptr[root1] > g_sample_ptr[root2])
			{
				//UpdateCenterOfMassAndEffectiveRadius(root2, root1);
				g_sample_ptr[root2] += g_sample_ptr[root1];
				g_sample_ptr[root1] = root2;
				UpdateBigClusterRoot(root2);
			}
			else
			{
				//UpdateCenterOfMassAndEffectiveRadius(root1, root2);
				g_sample_ptr[root1] += g_sample_ptr[root2];
				g_sample_ptr[root2] = root1;
				UpdateBigClusterRoot(root1);
			}
			if (big_clu_size[0] - size_of_old_biggest_clu > difference_of_biggest_clu)
				difference_of_biggest_clu = big_clu_size[0] - size_of_old_biggest_clu;
			if (kNumSample - big_clu_size[0] < difference_of_biggest_clu)
				break;
			/*
			distance_between_two_biggest_clu = 0;
			for (int j = 0; j < g_num_ising_site; j++) {
				distance_between_two_biggest_clu += (center_of_clu[big_clu_root[0]][j] -\
					center_of_clu[big_clu_root[1]][j])*\
					(center_of_clu[big_clu_root[0]][j] -\
					center_of_clu[big_clu_root[1]][j]);
			}
			distance_between_two_biggest_clu = sqrt(distance_between_two_biggest_clu);
			*/
		}
		/*
		if (big_clu_size[1] != 0) {
            outfile_clu_dis << big_clu_size[0] << " "  <<\
                effective_radius[big_clu_root[0]] << " ";
            outfile_clu_dis << big_clu_size[1] << " "  <<\
                effective_radius[big_clu_root[1]] << " " <<\
				distance_between_two_biggest_clu;
		}
		*/

		//fprintf(fp, "%lf %5d %5d %5d %5d %5d\n", (double)(i + 1) /kNumSample,\
		//	big_clu_size[0],big_clu_size[1], big_clu_size[2], big_clu_size[3], big_clu_size[4]);
	}
	//outfile_clu_dis.close();

	return difference_of_biggest_clu;
}

//Two schemes(functions) for calculating networks
//CalculateDelta: use InitializeOptimized
//CalculateDeltaMC: use InitialzeMC and MC networks 
int CalculateDelta()
{
	int difference_of_biggest_clu = 0;
	AllocateMemoryToArray();
	if (g_myrank == 0)
	{
		g_current_time = clock();
		cout << "memory finished-----process" << g_myrank << "----time:";
		cout << (double)(g_current_time - g_start_time) / CLOCKS_PER_SEC << endl;
	}
	InitializeOptimized();
	if (g_myrank == 0)
	{
		g_current_time = clock();
		cout << "initialization finished-----process" << g_myrank << "----time:";
		cout << (double)(g_current_time - g_start_time) / CLOCKS_PER_SEC << endl;
	}
	QuicksortDist(distances_sample, 0, num_add_bonds - 1);
	if (g_myrank == 0)
	{
		g_current_time = clock();
		cout << "sort finished-----process" << g_myrank << "----time:";
		cout << (double)(g_current_time - g_start_time) / CLOCKS_PER_SEC << endl;
	}
	//PrintDistance();
	difference_of_biggest_clu += BuildNetwork();
	if (g_myrank == 0)
	{
		g_current_time = clock();
		cout << "build finished-----process" << g_myrank << "----time:";
		cout << (double)(g_current_time - g_start_time) / CLOCKS_PER_SEC << endl;
	}
	DeleteMemoryToArray();
	return difference_of_biggest_clu;
}

double CalculateDeltaMC()
{
	double difference_of_biggest_clu = 0.0;
	int repeat_time = 20;
	for (int i = 0; i < repeat_time; i++)
	{
		AllocateMemoryToArrayMC();
		InitializeMC();
		QuicksortDist(distances_sample, 0, num_add_bonds - 1);
		//PrintDistance();
		difference_of_biggest_clu += BuildNetwork();
		DeleteMemoryToArray();
		if (g_myrank == 0)
		{
			g_current_time = clock();
			cout << i << "/" << repeat_time << endl;
			cout << "network finished-----process" << g_myrank << "----time:";
			cout << (double)(g_current_time - g_start_time) / CLOCKS_PER_SEC << endl;
		}
	}
	return difference_of_biggest_clu / repeat_time;
}