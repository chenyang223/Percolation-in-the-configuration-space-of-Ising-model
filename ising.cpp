#include "head.h"

int sample_count = 0;
double T_c;  //critical temperature
double *prob;   //probability for metropolis
double p_wolff; // probability for wolff

int **nearest_neighbor; // array to store the index of nearest neighbor

//calculate the nearest_neighbor for different g_ising_mode
void ComputeNeighbors()
{
    nearest_neighbor = new int *[g_num_ising_site];
    for (int i = 0; i < g_num_ising_site; i++)
        nearest_neighbor[i] = new int[g_num_nearest_neighbor];
    
    switch (g_ising_mode)
    {
    case 4:
    {
        for (int i = 0; i < g_num_ising_site; i++)
        {
            nearest_neighbor[i][0] = (i + 1) % g_num_ising_site;                             //right
            nearest_neighbor[i][1] = (i + g_num_ising_site - 1) % g_num_ising_site;            //left
            nearest_neighbor[i][2] = (i + g_length_ising) % g_num_ising_site;                  //down
            nearest_neighbor[i][3] = (i + g_num_ising_site - g_length_ising) % g_num_ising_site; //up
            if (i % g_length_ising == 0)
                nearest_neighbor[i][1] = i + g_length_ising - 1;
            if ((i + 1) % g_length_ising == 0)
                nearest_neighbor[i][0] = i - g_length_ising + 1;
        }
        break;
    }
    case 3:
    {
        for (int i = 0; i < g_num_ising_site; i++)
        {
            int row = i / g_length_ising;
            int column = i % g_length_ising;
            nearest_neighbor[i][0] = (column + 1) % g_length_ising + row * g_length_ising;
            nearest_neighbor[i][1] = (column - 1 + g_length_ising) % g_length_ising +
                                     row * g_length_ising;
            if ((column - row + g_length_ising) % 2 == 0)
                nearest_neighbor[i][2] = (i + g_num_ising_site - g_length_ising) % g_num_ising_site;
            else
                nearest_neighbor[i][2] = (i + g_length_ising) % g_num_ising_site;
        }
        break;
    }
    case 6:
    {
        for (int i = 0; i < g_num_ising_site; i++)
        {
            int row = i / g_length_ising;
            int column = i % g_length_ising;
            //same row
            nearest_neighbor[i][0] = (column + 1) % g_length_ising + row * g_length_ising;
            nearest_neighbor[i][1] = (column - 1 + g_length_ising) % g_length_ising +
                                     row * g_length_ising;
            //upper row
            nearest_neighbor[i][2] = column +
                                     (row - 1 + g_length_ising) % g_length_ising * g_length_ising;
            nearest_neighbor[i][3] = (column + 1) % g_length_ising +
                                     (row - 1 + g_length_ising) % g_length_ising * g_length_ising;
            //lower row
            nearest_neighbor[i][4] = column +
                                     (row + 1) % g_length_ising * g_length_ising;
            nearest_neighbor[i][5] = (column - 1 + g_length_ising) % g_length_ising +
                                     (row + 1) % g_length_ising * g_length_ising;
        }
        break;
    }
    case -6:
    {
        int square = g_length_ising * g_length_ising;
        for (int i = 0; i < g_num_ising_site; i++)
        {
            int height = i / (square);
            int row = (i % (square)) / g_length_ising;
            int column = (i % (square)) % g_length_ising;
            //same plane
            nearest_neighbor[i][0] = (column + 1) % g_length_ising + row * g_length_ising + height * square;
            nearest_neighbor[i][1] = (column - 1 + g_length_ising) % g_length_ising +
                                     row * g_length_ising + height * square;
            nearest_neighbor[i][2] = column + height * square +
                                     (row - 1 + g_length_ising) % g_length_ising * g_length_ising;
            nearest_neighbor[i][3] = column + height * square +
                                     (row + 1) % g_length_ising * g_length_ising;
            //different plane
            nearest_neighbor[i][4] = column + row * g_length_ising +
                                     (height + 1) % g_length_ising * square;
            nearest_neighbor[i][5] = column + row * g_length_ising +
                                     (height - 1 + g_length_ising) % g_length_ising * square;
        }
		break;
    }
    default:
    {
        cout << "Unavailable mode of Ising model." << endl;
        break;
    }
    }
}

void Wolff()
{
    int *stack = new int[g_num_ising_site]; //LIFO
    int stack_pointer = 1;
    int oldspin, newspin;

    int current;
    int seed = RandInt(g_num_ising_site);
    oldspin = g_ising_sites[seed];
    g_ising_sites[seed] *= -1;
    newspin = g_ising_sites[seed];
    stack[0] = seed;

    while (stack_pointer)
    {
        stack_pointer--;
        current = stack[stack_pointer];

        for (int i = 0; i < g_num_nearest_neighbor; i++)
        {
            if (g_ising_sites[nearest_neighbor[current][i]] == oldspin)
                if (RandDouble(1) < p_wolff)
                {
                    stack[stack_pointer] = nearest_neighbor[current][i];
                    g_ising_sites[stack[stack_pointer]] *= -1;
                    stack_pointer++;
                }
        }
    }
    delete[] stack;
}

void Metropolis()
{
    int flip = RandInt(g_num_ising_site);
    g_ising_sites[flip] *= -1;
    int spin_sum_nn = 0; //n.n. spin sum
    for (int i = 0; i < g_num_nearest_neighbor; i++)
        spin_sum_nn += g_ising_sites[nearest_neighbor[flip][i]];
    spin_sum_nn *= -g_ising_sites[flip];
    if (spin_sum_nn > 0 && prob[spin_sum_nn] < RandDouble(1))
        g_ising_sites[flip] *= -1; //no flipping
}

void InitializeIsing()
{
    switch (g_ising_mode)
    {
    case 4:
        T_c = 2.269;
        break;
    case 6:
        T_c = 3.643;
        break;
    case 3:
        T_c = 1.518;
        break;
    case -6:
        T_c = 4.51;
        break;
    default:
        cout << "unavaliable ising mode" << endl;
        break;
    }
    prob = new double[g_num_nearest_neighbor + 1];
    for (int i = 0; i < g_num_ising_site; i++)
        g_ising_sites[i] = (rand() % 2) ? 1 : -1;
    //g_ising_sites[i] = 1;
    for (int i = 0; i < g_num_nearest_neighbor + 1; i++)
        prob[i] = exp(-2 * g_J * i / g_T); //Metropolis possible value
    p_wolff = 1 - exp(-2 * g_J / g_T);     //wolff probability
}

void DeleteMemory()
{
    delete[] prob;
    for (int i = 0; i < g_num_nearest_neighbor; i++)
        delete[] nearest_neighbor[i];
    delete[] nearest_neighbor;
}

void SampleConfiguration()
{
    InitializeIsing();
    ComputeNeighbors();
    sample_count = 0; 

    int steps1; //steps for reaching equilibrium
    int steps2; //steps after equilibrium
    int sample_gap;

    if (g_T > T_c - 0.2 && g_T < T_c + 0.2)
    {
        if (g_T > T_c)
        {
            steps1 = 20 * (int)(2 * sqrt(g_length_ising) * (200 * abs(g_T - T_c) + 1));
            sample_gap = (int)(2 * sqrt(g_length_ising) * (200 * abs(g_T - T_c) + 1));
        }
        else
        {
            steps1 = 20 * (int)(2 * sqrt(g_length_ising) * (20 * abs(g_T - T_c) + 1));
            sample_gap = (int)(2 * sqrt(g_length_ising) * (20 * abs(g_T - T_c) + 1));
        }
        steps2 = sample_gap * (kNumSample - 1) + 1;

        for (int i = 0; i < steps1; i++)
            Wolff();
        for (int i = 0; i < steps2; i++)
        {
            Wolff();
            if (i % sample_gap == 0)
            {
                for (int j = 0; j < g_num_ising_site; j++)
                    g_sample[sample_count][j] = g_ising_sites[j];
                sample_count++;
                //cout << sample_count << endl;
            }
        }
    }
    else
    {
        steps1 = 10 * g_length_ising * g_num_ising_site;
        sample_gap = g_num_ising_site * g_length_ising;
        steps2 = sample_gap * (kNumSample - 1) + 1;

        for (int i = 0; i < steps1; i++)
            Metropolis();
        for (int i = 0; i < steps2; i++)
        {
            Metropolis();
            if (i % sample_gap == 0)
            {
                for (int j = 0; j < g_num_ising_site; j++)
                    g_sample[sample_count][j] = g_ising_sites[j];
                sample_count++;
            }
        }
    }
    DeleteMemory();
}

//monte carlo for test
void UpdateHm(double &H, double &m)
{
    H = 0;
    m = 0;
    for (int i = 0; i < g_num_ising_site; i++)
    {
        m += g_ising_sites[i];
        H += -g_J * g_ising_sites[i] * (g_ising_sites[nearest_neighbor[i][0]] + g_ising_sites[nearest_neighbor[i][2]]);
    }
}

//available for Monte Carlo simulation
void MonteCarloWolff()
{
    ComputeNeighbors();

    int steps1; //to reach equilibrium
    int steps2; //after reaching equilibrium
    int sample_gap;

    double H, m;

    double T_min = 4.45;
    double delta_T = 0.01;
    int const T_num = 11;
    double Tmax = T_min + T_num * delta_T;
    double Tcurrent = T_min;

    double H_avg[T_num], H_square_avg[T_num];
    double C[T_num];                                           //heat capacity
    double m_avg[T_num], m_sqare_avg[T_num], m_4th_avg[T_num]; //total spin
    double Chi[T_num];                                         //susceptability

    string file_name;
    stringstream convert_stream;
    convert_stream << g_length_ising << "_" << g_ising_mode << "_";
    convert_stream >> file_name;
    file_name.append("L T m c Chi Binder.txt");

    ofstream outfile;
    outfile.open(file_name.c_str());

    for (int i = 0; i < T_num; Tcurrent += delta_T, i++)
    {
        g_T = Tcurrent;
        InitializeIsing();

        if (g_T > T_c)
        {
            steps1 = 30 * (int)(2 * sqrt(g_length_ising) * (20 * abs(g_T - T_c) + 1));
            sample_gap = (int)(2 * sqrt(g_length_ising) * (20 * abs(g_T - T_c) + 1));
        }
        else
        {
            steps1 = 30 * (int)(2 * sqrt(g_length_ising)) * (20 * abs(g_T - T_c) + 1);
            sample_gap = (int)(2 * sqrt(g_length_ising)) * (20 * abs(g_T - T_c) + 1);
        }
        steps2 = sample_gap * (kNumSample - 1) + 1;

        H_avg[i] = 0;
        H_square_avg[i] = 0;
        m_avg[i] = 0;
        m_sqare_avg[i] = 0;
        m_4th_avg[i] = 0;

        for (int j = 0; j < steps1; j++)
            Wolff(); //reaching equilibrium
        for (int j = 0; j < steps2; j++)
        {
            Wolff();
            if (j % sample_gap == 0)
            {
                UpdateHm(H, m);
                m /= g_num_ising_site;
                H_avg[i] += H;
                H_square_avg[i] += H * H;
                m_avg[i] += fabs(m);
                m_sqare_avg[i] += m * m;
                m_4th_avg[i] += m * m * m * m;
            }
        }

        H_avg[i] /= (double)kNumSample;
        m_avg[i] /= (double)kNumSample;
        H_square_avg[i] /= (double)kNumSample;
        m_sqare_avg[i] /= (double)kNumSample;
        m_4th_avg[i] /= (double)kNumSample;

        C[i] = (H_square_avg[i] - H_avg[i] * H_avg[i]) / g_T / g_T / (g_length_ising * g_length_ising);
        Chi[i] = (m_sqare_avg[i] - m_avg[i] * m_avg[i]) / g_T * g_length_ising * g_length_ising;

        cout << g_T << endl;
        outfile << g_length_ising << " " << g_T << " " << fabs(m_avg[i]) << " ";
        outfile << C[i] << " " << Chi[i] << " ";
        outfile << m_4th_avg[i] / (m_sqare_avg[i] * m_sqare_avg[i]) << endl;
    }
    DeleteMemory();
}
