#include <iostream>
#include <fstream>
#include <random>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <list>
using namespace std;


/*
We will be using the vector for the location
*/
#include <stdio.h>
#include <vector>

typedef vector<int> int_vector;


const double TL = 1.0;
const double TR = 2.0;


struct interaction
{
    double time;
    vector<int> location;
    interaction* left;
    interaction* right;
};

void print_v(double* Array, int size)
{
    for(int i = 0; i < size; i++)
        cout<< Array[i] << "  ";
    cout<<endl;
}

int_vector find_min(interaction* Link)
{
        double tmp = Link->time;
        
        /*
        Initializing the starting point
        */
        
        interaction* pt = Link;
        interaction* tmp_pt = Link;
        
        /*
        Keep moving right until it reaches TR
        */
        while(tmp_pt != NULL)
        {
        	/*
        	Check if tmp is smaller than tmp_pt
        	*/
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
    /*
    Remember that interaction** represents the buckets
    */
	
	/*
	Check if the array in a bucket is empty or the point is empty
	*/
	if( *Link == NULL || pt == NULL)
        return;
        
    /*
    Check if the bucket has only one interaction
    */
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
    
    /*
    Now, isolate the pt that we just extracted
    */
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
        cout<<" location: "<< tmp->location <<" time: " << tmp->time << "  " ;
        tmp = tmp->right;
    }
    cout<<endl;
}


void big_step_distribute(interaction** &clock_time_in_step, interaction* time_array, 
	const int N1, const int N2, const double small_tau, const int ratio, const int Step)
//distribute clock times of a big step into vectors that represent small steps. 
//If the clock time is bigger than a big tau, then it is arranged in the right location
{
    for(int i = 0; i < N1; i++)
    {
    	for(int j = 0; j < N2; j++)
    	{
    		
        	int tmp;
        	if(time_array[i][j].time > (Step + 1)*ratio*small_tau)
        	{
            	tmp = ratio;
        	}
        	else
        	{
            	tmp = int((time_array[i][j].time - ratio*small_tau*Step)/small_tau);
        	}
        
        	push_front(&clock_time_in_step[tmp], &time_array[i][j]);
		}
        
    }
}


/*
void big_step_distribute(interaction** &clock_time_in_step, interaction* time_array, 
	const int N, const double small_tau, const int ratio, const int Step)
//distribute clock times of a big step into vectors that represent small steps. 
//If the clock time is bigger than a big tau, then it is arranged in the right location
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
*/


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
    }
    //    cout<<"start to move "<< pt->location << " from " << old_level << " to " << new_level<<endl;
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



void update(interaction** &clock_time_in_step, const int level, const int N, const int N2, const double small_tau,
 const int ratio, const int Step, interaction* time_array, double *energy_array, 
 uniform_real_distribution<double> &u, mt19937 &mt, int &count)
//update clock_time_in_step[level]
{
    double next_time = (Step*ratio + level + 1)*small_tau;
    
    int_vector min_loc = find_min(clock_time_in_step[level]);
    
    min_loc_row = min_loc.front();
    min_loc_col = min_loc.back();
    
    interaction* pt = &time_array[min_loc_col][min_loc_row];
    
    double current_time = pt->time;
    
//    cout<<"at level " << level << endl;

    while(current_time < next_time)
    {
        count++;

        //Step 1: update min interaction and energy
        double total_energy = energy_array[min_loc_col][min_loc_row] + energy_array[min_loc_col + 1][min_loc_row] 
			+ energy_array[min_loc_col][min_loc_row + 1] + energy_array[min_loc_col][min_loc_row - 1];
        
        double tmp_double = -log(u(mt))/sqrt(total_energy);
        
//        cout<<"added time = " << tmp_double << endl;

        double old_e_left = energy_array[min_loc_col][min_loc_row];
        double old_e_right = energy_array[min_loc_col + 1][min_loc_row];
        double old_e_up = energy_array[min_loc_col][min_loc_row + 1];
        double old_e_down = energy_array[min_loc_col][min_loc_row - 1];
        
        double tmp_rnd_uni = u(mt);
        
        if(min_loc_col == 0)
        {
        	/*if(min_loc_row == 1 || min_loc_row == N2)
        	{
        		total_energy = old_e_right - log(u(mt))*TL;
			}*/
			
			total_energy = old_e_right - log(u(mt))*TL;
			
			/*
			else
			{
				total_energy = old_e_right + old_e_down + old_e_up - log(u(mt))*TL;
			}*/
            // total_energy = old_e_right - log(u(mt))*TL;
        }
        
        if(min_loc_col == N)
        {
            if(min_loc_row == 1)
        	{
        		total_energy = old_e_left + old_e_up - log(u(mt))*TL;
			}
        	
        	else if(min_loc_row == N2)
        	{
        		total_energy = old_e_left + old_e_down - log(u(mt))*TL;
			}
			
			else
			{
				total_energy = old_e_left + old_e_down + old_e_up - log(u(mt))*TL;
			}
			// total_energy = old_e_left  -log(u(mt))*TR;
        }
        
        
        if(min_loc_col != 0)
        {
            energy_array[min_loc_col][min_loc_row] = tmp_rnd_uni*total_energy;
//            E_avg[min_loc - 1] += (current_time - last_update[min_loc - 1])*old_e_left;
//            last_update[min_loc - 1] = current_time;
        }
        
        if(min_loc_col != N)
        {
            energy_array[min_loc_col + 1][min_loc_row] = (1 - tmp_rnd_uni)*total_energy;//update energy
//            E_avg[min_loc] += (current_time - last_update[min_loc])*old_e_right;
//            last_update[min_loc] = current_time;
        }
        
        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step,  current_time + tmp_double);
        
        
        //Step 2: update other interactions
        if(min_loc != 0)
        {
            pt = &time_array[min_loc - 1];
            tmp_double = (pt->time - current_time)*sqrt(energy_array[min_loc - 1] 
            	+ old_e_left)/sqrt(energy_array[min_loc - 1] + energy_array[min_loc]) + current_time;
            move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
        }
        
        if(min_loc != N)
        {
            pt = &time_array[min_loc + 1];
            tmp_double = (pt->time - current_time)*sqrt(energy_array[min_loc + 2] 
            	+ old_e_right)/sqrt(energy_array[min_loc + 2] + energy_array[min_loc + 1]) + current_time;
            move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
            
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
        }
        
        
    }
    
}

int main(int argc, char** argv)
{
    struct timeval t1, t2;
    ofstream myfile;
    myfile.open("HL_KMP.txt", ios_base::app);

    int N = 10;
    int N2 = 10;
    
    if(argc > 1)
    {
        N = strtol(argv[1], NULL,10 );
        N2 = strtol(argv[1], NULL,10 );
    }
    
    double big_tau = 0.2; //big time step of tau leaping
    const int ratio = int(N/10); //ratio of big step and small step
    double small_tau = big_tau/double(ratio); //small time step
    
    double* energy_array = new double[N + 2][N2 + 2];
    double *E_avg = new double[N][N2];
    double* last_update = new double[N][N2];
    
    for(int i = 0; i <N; i++)
    {
        for(int j = 0; j < N2; j++)
        {
        	// E_avg[i][j] = 0;
        	last_update[i][j] = 0;
		}
    }
    
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> u(0,1);
    
    for(int k = 0; k < N2; k++)
    {
    	energy_array[0][k] = TL;
    	energy_array[N+1][k] = TR;
	}
	
	for(int j = 0; j < N; j++)
	{
		energy_array[j][0] = NULL;
		energy_array[j][N2 + 1] = NULL;
	}

    
    for(int n = 1; n < N+1; n++)
    {
    	for(int m = 1; m < N2 + 1; m++)
    	{
    		energy_array[n][m] = 1;
		}
	}
        
        
    interaction* time_array = new interaction[N + 1][N2 + 1];
    
    for(int n = 0; n < N+1; n++)
    {
    	for(int m = 0; m < N2 + 1; m++)
    	{
    		time_array[n][m].time = -log(u(mt))/sqrt(energy_array[n][m] + energy_array[n + 1][m]);
        	time_array[n][m].location = {n, m};
        	time_array[n][m].left = NULL;
        	time_array[n][m].right = NULL;
		}
    }
    
    int count = 0;

    interaction** clock_time_in_step = new interaction*[ratio + 1];//each element in the array is the head of a list
    for(int i = 0; i < ratio + 1; i++)
    {
        clock_time_in_step[i] = NULL;
    }
    
    gettimeofday(&t1,NULL);
    
//    big_step_distribute(clock_time_in_step,time_array,N+1,small_tau,ratio,0);

    
    
    int Step = 50;
   
    for(int out_n = 0; out_n < Step; out_n++)
    {
//        cout<<"at Step "<<out_n<<endl;
        /*
        for(int i = 0; i <= ratio; i++)
        {
            cout<<"at level " << i << "   ";
            print_list( clock_time_in_step[i]);
        }
         */
        big_step_distribute(clock_time_in_step,time_array,N+1,small_tau,ratio,out_n);
        
        for(int in_n = 0; in_n < ratio; in_n++)
        {
            
            if(clock_time_in_step[in_n]!= NULL)
            {
                update(clock_time_in_step, in_n, N, small_tau, ratio, 
                	out_n, time_array, energy_array, u, mt, count);
            }
//            print_v(energy_array,N+2);
        }
        
        /*
        for(int i = 0; i <= ratio; i++)
        {
            cout<<"at level " << i << "   ";
            print_list( clock_time_in_step[i]);
        }
        */
        
        clock_time_in_step[ratio] = NULL;
        
        
        

        
    }
    
    /*
    for(int i = 0; i < N; i++)
    {
        E_avg[i]= E_avg[i]/(Step*big_tau);
    }
    
    print_v(E_avg, N);
    */
    gettimeofday(&t2, NULL);
    delete[][] energy_array;
    delete[][] E_avg;
    delete[][] time_array;
    delete[] clock_time_in_step;

    double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
                    t2.tv_usec - t1.tv_usec) / 1.e6;
    
	cout << "total CPU time = " << delta <<endl;
    cout<<" N = "<<N <<endl;
    cout<<" N2 = "<<N2 <<endl;
    
    //cout<<"seconds per million event is "<< 1000000*delta/double(count)<<endl;
    //myfile<<" N = "<<N <<endl;
    
    myfile<< 1000000*delta/double(count)<<"  ";
    myfile.close();

 
 
    
}


