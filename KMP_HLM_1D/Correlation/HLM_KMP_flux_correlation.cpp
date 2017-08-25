#include <iostream>
#include <fstream>
#include <random>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <list>
//#include <trng/yarn2.hpp>
//#include <trng/uniform01_dist.hpp>
using namespace std;

typedef list<double> d_List;

const double TL = 1.0;
const double TR = 10.0;

enum Rate_Func {
    rf1,
    rf2
};

struct interaction
{
    double time;
    int location;
    int level;
    interaction* left;
    interaction* right;
};

struct correlation_list
{
    double previous_time;
    double energy_sum;
};

Rate_Func rate_func = rf1;

double rate_function(double x, double y) 
{
    if(x >= 0 && y >= 0)
        return (rate_func == rf1) ? sqrt(x*y/(x+y)) : sqrt(x+y);
    else {
        cout<<"error, sqrt of negative number! "<<endl;
        return 0;
    }
}

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
        cout<<" location: "<< tmp->location <<" time: " << tmp->time << "  " ;
        tmp = tmp->right;
    }
    cout<<endl;
}


void big_step_distribute(interaction** &clock_time_in_step, interaction* time_array, const int N, 
    const double small_tau, const int ratio, const int Step)
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



void update(interaction** &clock_time_in_step, const int level, const int N, const double small_tau, 
    const int ratio, const int Step, interaction* time_array, double *energy_array, 
    uniform_real_distribution<double> &u, mt19937 &mt, int &count, correlation_list &left_cell, 
    correlation_list &right_cell, correlation_list &LR_multiplied)
//update clock_time_in_step[level]
{
    double next_time = (Step*ratio + level + 1)*small_tau;
    int min_loc = find_min(clock_time_in_step[level]);
    interaction* pt = &time_array[min_loc];
    double current_time = pt->time;

    double time_multiple1, time_multiple2, time_multiple3;
    int mid_point = int(N/2);

    while(current_time < next_time)
    {
        count++;

        //Step 1: update min interaction and energy
        double total_energy = energy_array[min_loc] + energy_array[min_loc + 1];
        double tmp_double = - log(1 - u(mt))/sqrt(total_energy);
//        cout<<"added time = " << tmp_double << endl;
        double old_e_left = energy_array[min_loc];
        double old_e_right = energy_array[min_loc + 1];
        double tmp_rnd_uni = u(mt);

        //int mid_point = int(N/2);

        if(min_loc == 0)
        {
            total_energy = old_e_right - log(1 - u(mt))*TL;
        }
        if(min_loc == N)
        {
            total_energy = old_e_left - log(1 - u(mt))*TR;
        }

        if(min_loc != 0)
        {
            energy_array[min_loc] = tmp_rnd_uni*total_energy;
        }
        if(min_loc != N)
        {
            energy_array[min_loc + 1] = (1 - tmp_rnd_uni)*total_energy;
        }


        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step,  current_time + tmp_double);
        
        //double time_multiple1, time_multiple2, time_multiple3;

        // left and right
        if(min_loc == mid_point)
        {
            //cout << "1" << endl;
            time_multiple1 = current_time - left_cell.previous_time;
            time_multiple2 = current_time - right_cell.previous_time;
            time_multiple3 = current_time - LR_multiplied.previous_time;
            // cout << time_multiple1 << " " << time_multiple2 << " " << time_multiple3 << endl;
            left_cell.previous_time = current_time;
            right_cell.previous_time = current_time;
            LR_multiplied.previous_time = current_time;

            left_cell.energy_sum += time_multiple1 * energy_array[mid_point];
            right_cell.energy_sum += time_multiple2 * energy_array[mid_point + 1];
            LR_multiplied.energy_sum += time_multiple3 * energy_array[mid_point] * energy_array[mid_point + 1];
            // left_cell.energy_sum += time_multiple1*old_e_left;
            // right_cell.energy_sum += time_multiple2*old_e_right;
            // LR_multiplied.energy_sum += time_multiple3*old_e_left*old_e_right;
        }

        // right
        else if(min_loc == mid_point + 1)
        {
            //cout << "2" << endl;
            time_multiple2 = current_time - right_cell.previous_time;
            time_multiple3 = current_time - LR_multiplied.previous_time;
            // cout << time_multiple2 << " " << time_multiple3 << endl;
            right_cell.previous_time = current_time;
            LR_multiplied.previous_time = current_time;

            right_cell.energy_sum += time_multiple2 * energy_array[mid_point + 1];
            LR_multiplied.energy_sum += time_multiple3 * energy_array[mid_point] * energy_array[mid_point + 1];
            // right_cell.energy_sum += time_multiple2*old_e_right;
            // LR_multiplied.energy_sum += time_multiple3*old_e_left*old_e_right;
        }

        // left
        else if(min_loc == mid_point - 1)
        {
            //cout << "3" << endl;
            time_multiple1 = current_time - left_cell.previous_time;
            time_multiple3 = current_time - LR_multiplied.previous_time;
            // cout << time_multiple1 << " " << time_multiple3 << endl;
            left_cell.previous_time = current_time;
            LR_multiplied.previous_time = current_time;

            left_cell.energy_sum += time_multiple1 * energy_array[mid_point];
            LR_multiplied.energy_sum += time_multiple3 * energy_array[mid_point] * energy_array[mid_point + 1];
            // left_cell.energy_sum += time_multiple1*old_e_left;
            // LR_multiplied.energy_sum += time_multiple3*old_e_left*old_e_right;
        }

        //cout << "left: " << left_cell.energy_sum << " right: " << right_cell.energy_sum << " LR: " << LR_multiplied.energy_sum  << endl;


        //Step 2: update other interactions
        if(min_loc != 0)
        {
            pt = &time_array[min_loc - 1];
            tmp_double = (pt->time - current_time)*sqrt(energy_array[min_loc - 1] + old_e_left)/
                sqrt(energy_array[min_loc - 1] + energy_array[min_loc]) + current_time;
            move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
        }
        if(min_loc != N)
        {
            pt = &time_array[min_loc + 1];
            tmp_double = (pt->time - current_time)*sqrt(energy_array[min_loc + 2] + old_e_right)/sqrt(energy_array[min_loc + 2] + energy_array[min_loc + 1]) + current_time;
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
            current_time = next_time + 100;
        }
        
        
    }
    
}

int main(int argc, char** argv)
{
    struct timeval t1, t2;
    ofstream myfile;
    myfile.open("HL_KMP.txt", ios_base::app);

    int N = 100;

    if(argc > 1)
    {
        N = strtol(argv[1], NULL,10 );
    }

    double big_tau = 0.2;//big time step of tau leaping
    const int ratio = int(N/10);//ratio of big step and small step
    double small_tau = big_tau/double(ratio);//small time step
    
    double* energy_array = new double[N+2];
    double *E_avg = new double[N];
    double* last_update = new double[N];

    for(int i = 0; i <N; i++)
    {
        E_avg[i] = 0;
        last_update[i] = 0;
    }

    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> u(0,1);

    energy_array[0] = TL;
    energy_array[N+1] = TR;

    for(int n = 1; n < N+1; n++)
    {
        energy_array[n] = 1;
    }

    interaction* time_array = new interaction[N+1];
    for(int n = 0; n < N+1; n++)
    {
        time_array[n].time = -log(u(mt))/sqrt(energy_array[n] + energy_array[n+1]);
        time_array[n].location = n;
        time_array[n].left = NULL;
        time_array[n].right = NULL;
    }

    int count = 0;

    interaction** clock_time_in_step = new interaction*[ratio + 1];
    //each element in the array is the head of a list
    for(int i = 0; i < ratio + 1; i++)
    {
        clock_time_in_step[i] = NULL;
    }
    
    gettimeofday(&t1,NULL);
    
    correlation_list X;
    X.energy_sum = 0;
    X.previous_time = 0;

    correlation_list Y;
    Y.energy_sum = 0;
    Y.previous_time = 0;

    correlation_list X_Y;
    X_Y.energy_sum = 0;
    X_Y.previous_time = 0;
    
    
    int Step = 100000;
   
    for(int out_n = 0; out_n < Step; out_n++)
    {
        big_step_distribute(clock_time_in_step, time_array, N+1 , small_tau, ratio, out_n);
        
        for(int in_n = 0; in_n < ratio; in_n++)
        {
            
            if(clock_time_in_step[in_n]!= NULL)
            {
                update(clock_time_in_step, in_n, N, small_tau, ratio, out_n, time_array, energy_array, u, 
                    mt, count, X, Y, X_Y);
            }
//            print_v(energy_array,N+2);
        }
        
        clock_time_in_step[ratio] = NULL;
        
    }

    cout << "E[X] = " << X.energy_sum/(Step*big_tau) << endl;
    cout << "E[Y] = " << Y.energy_sum/(Step*big_tau) << endl;
    cout << "E[XY] = " << X_Y.energy_sum/(Step*big_tau) << endl;
    cout << "Hence, Cov(X, Y) = " << (X_Y.energy_sum)/(Step*big_tau) - (X.energy_sum * Y.energy_sum)/(Step*big_tau*Step*big_tau) << endl;

    gettimeofday(&t2, NULL);
    delete[] energy_array;
    delete[] E_avg;
    delete[] time_array;
    delete[] clock_time_in_step;

    double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
                    t2.tv_usec - t1.tv_usec) / 1.e6;
    
//    cout << "total CPU time = " << delta <<endl;
    cout<<" N = "<<N <<endl;
    //cout<<"seconds per million event is "<< 1000000*delta/double(count)<<endl;
    //myfile<<" N = "<<N <<endl;
    myfile<< 1000000*delta/double(count)<<"  ";
    myfile.close();


}