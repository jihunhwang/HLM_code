/*
Now please do the following tests:

1, Mean energy flux for changing length of the chain. Let N change from 2 to 100.

2, Mean energy flux for changing difference of boundary temperatures. Fix TL = 1 and change TR from 1.5 to 10.0.

3, Mean energy flux for changing temperatures. Change TL from 1 to 100 and let TR be 1.10 times TL.

Compute the mean energy flux for each setting for 10 times and calculate the standard deviation. The standard deviation should be small comparing to the value of mean energy flux. If the computation task is too much for a laptop, you can ssh to my desktop “io.math.umass.edu”. The desktop has 8 cores and 16 threads.
*/


/*
This is a test file for test #1: Mean energy flux for changing length.
*/


#include <iostream>
#include <fstream>
#include <random>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <list>
#include <omp.h>
//#include <array>
//#include <iterator>

#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

using namespace std;

/*
template<typename T, std::size_t N>
constexpr std::size_t size(T(&)[N])
{
	return N;
}*/


const double TL = 1.0;
const double TR = 2.0;

struct interaction
{
    double time;
    int location;
    interaction* left;
    interaction* right;
};


void print_v(double* Array, int size)
{
    for(int i = 0; i < size; i++)
        cout<< Array[i] << "  ";
    cout<<endl;
}


int find_min(interaction* Link)
{
        double tmp = Link->time;
        interaction* pt = Link;
        interaction* tmp_pt = Link;
        while(tmp_pt != NULL)
        {
            if(tmp_pt->time < tmp)
            {
                tmp = tmp_pt->time;
                pt = tmp_pt;
            }
            
            tmp_pt = tmp_pt->right;
        }
        return pt->location;
}


void remove(interaction** Link, interaction* pt)
//remove pt from link but keep pt for future use
{
    if( *Link == NULL || pt == NULL)
        return;
    else if( pt->left == NULL && pt->right == NULL)
    {
        *Link = NULL;
        return;
    }
    else
    {
        if((*Link) == pt )
            *Link = pt->right;
        if( pt->right != NULL)
            pt->right->left = pt->left;
        if( pt->left != NULL)
            pt->left->right = pt->right;
    }
    pt->left = NULL;
    pt->right = NULL;
}


void push_front(interaction** Link, interaction* pt)
//push the interaction pointed by pt into the front of the lise
{
    pt->left = NULL;
    if(*Link == NULL)
    {
        pt->right = NULL;
        *Link = pt;
    }
    else
    {
        pt->right = *Link;
        (*Link)->left = pt;
        *Link = pt;
    }
}


void print_list(interaction* Link)
{
    interaction* tmp = Link;
    while(tmp!= NULL)
    {
        cout << " location: " << tmp->location << " time: " << tmp->time << "  " ;
        tmp = tmp->right;
    }
    cout<<endl;
}


void big_step_distribute(interaction** &clock_time_in_step, interaction* time_array, 
     const int N, const double small_tau, const int ratio, const int Step)
//distribute clock times of a big step into vectors that represent small steps. If the clock time is bigger than a big tau, then it is arranged in the right location
{
    for(int i = 0; i < N; i++)
    {
        int tmp;
        if(time_array[i].time > (Step + 1)*ratio*small_tau)
        {
            tmp = ratio;
        }
        else
        {
            tmp = int((time_array[i].time - ratio*small_tau*Step)/small_tau);
        }
        
        push_front(&clock_time_in_step[tmp], &time_array[i] );
        
    }
}



void move_interaction(interaction** &clock_time_in_step, interaction* pt, const double small_tau, 
     const int ratio, const int Step, const double new_time)
//move the interaction pointed by *pt from old bucket to new bucket
//n_move1: move without relinking pointers
//n_move2: move with relinking pointers
{
    double old_time = pt->time;
    int old_level, new_level;
    pt->time = new_time;
    
    if (old_time > (Step + 1)*ratio*small_tau )
    {
        old_level = ratio;
    }
    
    else
    {
        old_level = int( (old_time - Step*ratio*small_tau)/small_tau );
    }
    
    
    if ( new_time > (Step + 1)*ratio*small_tau )
    {
        new_level = ratio;
    }
    
    else
    {
        new_level = int( (new_time - Step*ratio*small_tau)/small_tau );
        if( new_level < 0 || new_level >= ratio)
        {
            cout<<"Error!!"<<endl;
            cout<<"start to move "<< pt->location << " from " << old_level << " to " << new_level<<endl;
            cout<<"error time is "<<new_time<<endl;
            
        }
    }
    
    
    if(old_level == new_level )
    {
        pt->time = new_time;
    }
    
    else
    {
        remove(&(clock_time_in_step[old_level]), pt);
        push_front(&(clock_time_in_step[new_level]), pt);
    }
}


double rate_function(double x, double y)
{
    return sqrt(x * y/(x + y));
}

/*double standard_deviation(double* profilearray, int Nthread)
{
	double prof_arr_avg = 0.0;
	//int arr_len = end(profilearray) - begin(profilearray);
	int arr_len = sizeof(profilearray);
	//int arr_len = Nthread;

	//cout << arr_len << endl;
	//cout << profilearray << endl;

	for(int a = 0; a < arr_len; a++)
	{
		prof_arr_avg += profilearray[a];
		cout << "array: " << profilearray[a] << endl;
	}

	prof_arr_avg = prof_arr_avg/arr_len;

	double std_return = 0.0;

	for(int a = 0; a < arr_len; a++)
	{
		std_return += (profilearray[a] - prof_arr_avg) * (profilearray[a] - prof_arr_avg);
	}

	std_return = sqrt(std_return/arr_len);

	return std_return;
}*/

void update(interaction** &clock_time_in_step, const int level, const int N, const double small_tau, 
     const int ratio, const int Step, interaction* time_array, double *energy_array, 
     trng::uniform01_dist<> &u, trng::yarn2 &r, int &count, double &flux_J)
//update clock_time_in_step[level]
{
    double next_time = (Step*ratio + level + 1)*small_tau;
    int min_loc = find_min(clock_time_in_step[level]);
    interaction* pt = &time_array[min_loc];
    double current_time = pt->time;
    
    //while(current_time < next_time)
    while(clock_time_in_step[level] != NULL)
    {
        count++;

        //Step 1: update min interaction and energy
        double total_energy = energy_array[min_loc] + energy_array[min_loc + 1];
        double tmp_double;// = -log(1 - u(r))/sqrt(total_energy);
        double rnd;
        double tmp_rate;

        do
        {
            rnd = u(r);
        }while(rnd < 1e-16 || rnd > 1 - 1e-16);

        tmp_rate = rate_function(energy_array[min_loc], energy_array[min_loc + 1]);

        if(tmp_rate < 1e-16|| std::isnan(tmp_rate))
        {
           tmp_rate = 1e-16;
           cout<<"energy ="<< energy_array[min_loc]<<" "<<energy_array[min_loc + 1]<<endl;
        }

        tmp_double = - log(rnd)/tmp_rate;

        double old_e_left = energy_array[min_loc];
        double old_e_right = energy_array[min_loc + 1];
        double tmp_rnd_uni = u(r);
        
        do{
            tmp_rnd_uni = u(r);
        }while(tmp_rnd_uni < 1e-16 || tmp_rnd_uni > 1 - 1e-16);

        double previous_energy = energy_array[min_loc];

        if(min_loc == 0)
        {
            total_energy = old_e_right - log(1 - u(r))*TL;
            previous_energy = energy_array[min_loc + 1];
        }
        
        if(min_loc == N)
        {
            total_energy = old_e_left - log(1 - u(r))*TR;
        }
        
        if(min_loc != 0)
        {
            energy_array[min_loc] = tmp_rnd_uni*total_energy;
        }
        
        if(min_loc != N)
        {
            energy_array[min_loc + 1] = (1 - tmp_rnd_uni)*total_energy;
        }
        
        double current_energy = energy_array[min_loc];

        if(min_loc == 0)
        {
        	current_energy = energy_array[min_loc + 1];
        }

        flux_J = flux_J + ((1 - 2*int(min_loc == 0))*(current_energy - previous_energy));

        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step,  current_time + tmp_double);
        
        
        //Step 2: update other interactions
        if(min_loc != 0)
        {
            pt = &time_array[min_loc - 1];
            tmp_double = (pt->time - current_time)*rate_function(energy_array[min_loc - 1], old_e_left)/
                       rate_function(energy_array[min_loc - 1], energy_array[min_loc]) + current_time;
            move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
        }
        if(min_loc != N)
        {
            pt = &time_array[min_loc + 1];
            tmp_double = (pt->time - current_time)*rate_function(energy_array[min_loc + 2], old_e_right)/
                       rate_function(energy_array[min_loc + 2], energy_array[min_loc + 1]) + current_time;
            move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
        }
        
        
        
        //Step 3: update current time
        if(clock_time_in_step[level] != NULL)
        {
            min_loc = find_min(clock_time_in_step[level]);
            pt = &time_array[min_loc];
            current_time = pt->time;
        }
        else
        {
            current_time = next_time + 1;
            break;
        }
        
    }
}



int main(int argc, char** argv)
{
    struct timeval t1, t2;
    ofstream myfile;
    ofstream myjimmy;
    myfile.open("HL_KMP.txt", ios_base::app);
    myjimmy.open("KMP_Flux_Test1-1.txt", ios_base::trunc);

    const int Step = 1000000;
    int N_thread = 4;
    
	cout << " " << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "   Test result of change of mean energy flux" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "   Test #1 - 1" << endl;
    cout << "   Control variable: Length of the KMP Chain (N)" << endl;
    cout << "	Rate function: sqrt(x * y/(x + y))" << endl;
    cout << "   Step = " << Step << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << " " << endl;
    cout << " " << endl;

    myjimmy << " " << endl;
    myjimmy << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    myjimmy << "   Test result of change of mean energy flux" << endl;
    myjimmy << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    myjimmy << "   Test #1 - 1" << endl;
    myjimmy << "   Control variable: Length of the KMP Chain (N)" << endl;
    myjimmy << "   Rate function: sqrt(x * y/(x + y))" << endl;
    myjimmy << "   Step = " << Step << endl;
    myjimmy << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    myjimmy << " " << endl;
    myjimmy << " " << endl;

    int N = 100;

    if(argc > 1)
    {
        N = strtol(argv[1], NULL,10 );
    }

	// big time step of tau leaping
    double big_tau = 0.2; 

    // ratio of big step and small step
    int ratio;// = int(N / 10); 
    
    // small time step
    double small_tau; //= big_tau/double(ratio); 

    //int N_thread = 4;
    //const int Step = 10000000;

    // parallel_flux is an array that stores the energy computed from each threads
    double* parallel_flux = new double[N_thread];
    //double parallel_flux[N_thread];

    // initialize the parallel_flux
    for(int i = 0; i < N_thread; i++)
    {
    	parallel_flux[i] = 0.0;
    }

 	// Give a change to N (length of the chain)
    for(N = 60; N <= 60; N++)
    {
    	ratio = (N < 10) ? 2 : int(N / 10);
        small_tau = (N < 10) ? 0.1 : big_tau/double(ratio);

	    // We are going to repeat each N's ten times, and compute the average of ten trials.
	    double* average_stored = new double[10];

	    // initialization of average_stored
	    for(int j = 0; j < 10; j++)
	    {
	    	average_stored[j] = 0.0;
	    }

	    // for loop for repeating ten times
	    // repeating a same computation ten times will decrease the standard deviation by 90%
	    for(int cnt = 0; cnt < 10; cnt++)
		{
		    
		    cout << "-----------------------------------" << endl;
			cout << " N = " << N << endl;
		    cout << "Trial#: " << cnt + 1 << endl;
		    cout << "Flux (J) from each four cores: " << endl;

		    myjimmy << "-----------------------------------" << endl;
		    myjimmy << " N = " << N << endl;
		    myjimmy << "Trial#: " << cnt + 1 << endl;
		    myjimmy << "Flux (J) from each four cores: " << endl;

		    // Re-initialize the parallel_flux in the case it's not zero
		    // If it's not zero, it'll be ended up being accumulated with previous values
		    //for(int i = 0; i < N_thread; i++)
    		//{
    		//	parallel_flux[i] = 0.0;
    		//}


			// OpenMP starts here
		    #pragma omp parallel num_threads(N_thread)
		    {
		    	// Number of threads
		    	int rank = omp_get_thread_num();

		    	double* energy_array = new double[N + 2];
		    	double *E_avg = new double[N];
		    	double* last_update = new double[N];

		    	// Re-initialize the parallel_flux in the case it's not zero
		    	// If it's not zero, it'll be ended up being accumulated with previous values
		    	for(int i = 0; i < N_thread; i++)
    			{
    				parallel_flux[i] = 0.0;
    			}

		    	for(int i = 0; i < N; i++)
		    	{
		        	E_avg[i] = 0;
		        	last_update[i] = 0;
		    	}

		    	// TRNG supports OpenMP.
		    	// other RNGs don't usually work with parallel programming
		    	trng::yarn2 r;
		        trng::uniform01_dist<> u;
		        r.seed(time(NULL));
		        r.split(N_thread, rank);

		    	energy_array[0] = TL;
		    	energy_array[N + 1] = TR;

		    	for(int n = 1; n < N + 1; n++)
		    	{
		        	energy_array[n] = 1;
		    	}

		    	interaction* time_array = new interaction[N+1];

		    	for(int n = 0; n < N + 1; n++)
		    	{
		        	time_array[n].time = -log(1 - u(r))/rate_function(energy_array[n], energy_array[n + 1]);
		        	time_array[n].location = n;
		        	time_array[n].left = NULL;
		        	time_array[n].right = NULL;
		    	}

		    	//cout << "4" << endl;
		    	int count = 0;
		    
			    // each element in the array is the head of a list
			    interaction** clock_time_in_step = new interaction*[ratio + 1]; 

			    for(int i = 0; i < ratio + 1; i++)
			    {
			        clock_time_in_step[i] = NULL;
			    }
		    	//cout << "5" << endl;

			    gettimeofday(&t1,NULL);
			    
			    //int Step = 100000;
			    
			    //double flux_J = 0.0;
			
				//cout << "6" << endl;

			    for(int out_n = 0; out_n < Step; out_n++)
			    {
			        big_step_distribute(clock_time_in_step, time_array, N + 1, small_tau, ratio, out_n);
			        
			        for(int in_n = 0; in_n < ratio; in_n++)
			        {
			            if(clock_time_in_step[in_n]!= NULL)
			            {
			                update(clock_time_in_step, in_n, N, small_tau, ratio, out_n, time_array, energy_array, 
			                	u, r, count, parallel_flux[rank]);
			            }
			        }
			        clock_time_in_step[ratio] = NULL;
			    }
			    
			    
			    gettimeofday(&t2, NULL);
			    /*
			    delete[] energy_array;
			    delete[] E_avg;
			    delete[] time_array;
			    delete[] clock_time_in_step;*/

			    //double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u + t2.tv_usec - t1.tv_usec) / 1.e6;
			    //double large_T = 1000000*delta/double(count);

			    //cout << parallel_flux[rank]/large_T << endl;
			    cout << parallel_flux[rank] << endl;
		    }

	    	for(int i = 0; i < N_thread; i++)
	    	{
	    		myjimmy << parallel_flux[i] << endl;
	    	}

	    	double J_Avg = 0.0;
	    	double stand_dev = 0.0;

	    	for(int i = 0; i < N_thread; i++)
	    	{
	    		J_Avg += parallel_flux[i];
	    	}

			J_Avg = J_Avg/N_thread;

	    	cout << "Average: " << J_Avg << endl;
	    	myjimmy << "Average: " << J_Avg << endl;

	    	for(int i = 0; i < N_thread; i++)
	    	{
	    		stand_dev += (J_Avg - parallel_flux[i]) * (J_Avg - parallel_flux[i]);
	    	}

	    	//stand_dev = sqrt(stand_dev * 1/N_thread);
	    	cout << "Standard deviation: " << sqrt(stand_dev * 1/N_thread) << endl;
	    	myjimmy << "Standard deviation: " << sqrt(stand_dev * 1/N_thread) << endl;

            J_Avg = J_Avg/(Step*big_tau);
            cout << "Energy Flux: " << J_Avg << endl;
            myjimmy << "Energy Flux: " << J_Avg << endl;

	    	//cout << "Time" << large_T << endl;
	    	//myjimmy << "Time" << large_T << endl;

            //average_stored[cnt] = J_Avg/(N_thread*Step*big_tau);
	    	average_stored[cnt] = J_Avg;
	    	//cout << "J_Avg = " << J_Avg << endl;
	    }

	    // In this step, we are going to calculate the average and standard deviation
	    // of the energies we got from previous ten trials.
	    double aver_of_energy = 0.0;
	    double std_dev_of_ten = 0.0;


	    for(int b = 0; b < 10; b++)
	    {
	    	aver_of_energy += average_stored[b];
	    }

	    //aver_of_energy = aver_of_energy/(N_thread*Step*big_tau);
        aver_of_energy = aver_of_energy/10;
	    //cout << average_stored << endl;

	    cout << " " << endl;
	    cout << " " << endl;
	    cout << "************************ CONCLUSION *******************************" << endl;
	    cout << " " << endl;
	    cout << "When N = " << N << endl;
	    cout << " " << endl;
	    cout << "Energy Flux J = " << aver_of_energy << endl;

	    myjimmy << " " << endl;
	    myjimmy << " " << endl;
	    myjimmy << "************************ CONCLUSION *******************************" << endl;
	    myjimmy << " " << endl;
	    myjimmy << "When N = " << N << endl;
	    myjimmy << " " << endl;
	    myjimmy << "Energy Flux J = " << aver_of_energy << endl;


	    for(int c = 0; c < 10; c++)
	    {
	    	std_dev_of_ten += (aver_of_energy - average_stored[c]) * (aver_of_energy - average_stored[c]);
	    }

	    cout << "Standard Deviation = " << sqrt(std_dev_of_ten * 1/10) << endl;
	    cout << " " << endl;
	    cout << "*******************************************************************" << endl;
	    cout << " " << endl;
	    cout << " " << endl;

	    myjimmy << "Standard Deviation = " << sqrt(std_dev_of_ten * 1/10) << endl;
	    myjimmy << " " << endl;
	    myjimmy << "*******************************************************************" << endl;
	    myjimmy << " " << endl;
	    myjimmy << " " << endl;
	}

    //cout << "total CPU time = " << delta << endl;
    //cout << " N = " << N << endl;
    
    //cout<<"seconds per million event is "<< 1000000*delta/double(count)<<endl;
    //myfile << " N = " << N <<endl;
    //myfile << 1000000*delta/double(count) << "  ";
    myjimmy.close();
    myfile.close();
}

