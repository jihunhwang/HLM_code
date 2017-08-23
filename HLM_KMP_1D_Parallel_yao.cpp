#include <iostream>
#include <fstream>
#include <random>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <list>
#include <omp.h>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <array>
using namespace std;

const double TL = 1.0;
const double TR = 2.0;

struct interaction
{
    double time;
    int location;
    interaction* left;
    interaction* right;
};

struct thread_info {
    int current_rank = -1;
    double running_time = 0.0;
    int count = 0;
};

double rate_function(double x, double y) {
    
    return sqrt(x*y/(x+y));
    
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


void big_step_distribute(interaction** &clock_time_in_step, interaction* time_array, const int N, const double small_tau, const int ratio, const int Step)
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



void move_interaction(interaction** &clock_time_in_step, interaction* pt, const double small_tau, const int ratio, const int Step, const double new_time)
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



void update(interaction** &clock_time_in_step, const int level, const int N, const double small_tau, const int ratio, const int Step, interaction* time_array, double *energy_array, trng::uniform01_dist<> &u, trng::yarn2 &r, int &count, double &energy_integration)
//update clock_time_in_step[level]
{
    double next_time = (Step*ratio + level + 1)*small_tau;
    int min_loc = find_min(clock_time_in_step[level]);
    interaction* pt = &time_array[min_loc];
    double current_time = pt->time;
    double previous_energy = 0.0;
    double current_energy = 0.0;
//    cout<<"at level " << level << endl;
    while(current_time < next_time)
    {
        count++;

        //Step 1: update min interaction and energy
        double total_energy = energy_array[min_loc] + energy_array[min_loc + 1];
        double tmp_double = -log(1-u(r))/sqrt(total_energy);
//        cout<<"added time = " << tmp_double << endl;
        double old_e_left = energy_array[min_loc];
        double old_e_right = energy_array[min_loc + 1];
        double tmp_rnd_uni = u(r);

        previous_energy = energy_array[min_loc];
        if(min_loc == 0)
        {
            previous_energy = energy_array[min_loc + 1];
            total_energy = old_e_right  -log(1 - u(r))*TL;
        }
        if(min_loc == N)
        {
            total_energy = old_e_left  -log(1 - u(r))*TR;
        }
        if(min_loc != 0)
        {
            energy_array[min_loc] = tmp_rnd_uni*total_energy;
//            E_avg[min_loc - 1] += (current_time - last_update[min_loc - 1])*old_e_left;
//            last_update[min_loc - 1] = current_time;
        }
        if(min_loc != N)
        {
            energy_array[min_loc + 1] = (1 - tmp_rnd_uni)*total_energy;//update energy
//            E_avg[min_loc] += (current_time - last_update[min_loc])*old_e_right;
//            last_update[min_loc] = current_time;
        }
        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step,  current_time + tmp_double);
        current_energy = energy_array[min_loc];
        if(min_loc == 0) {
            current_energy = energy_array[min_loc + 1];
        }
        
        energy_integration += (1 - 2*int(min_loc == 0))*(current_energy - previous_energy);
        //Step 2: update other interactions
        if(min_loc != 0)
        {
            pt = &time_array[min_loc - 1];
            tmp_double = (pt->time - current_time)*rate_function(energy_array[min_loc - 1], old_e_left)/rate_function(energy_array[min_loc - 1], energy_array[min_loc]) + current_time;
            //tmp_double = (pt->time - current_time)*sqrt(energy_array[min_loc - 1] + old_e_left)/sqrt(energy_array[min_loc - 1] + energy_array[min_loc]) + current_time;
            move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
        }
        if(min_loc != N)
        {
            pt = &time_array[min_loc + 1];
            tmp_double = (pt->time - current_time)*rate_function(energy_array[min_loc + 2], old_e_right)/rate_function(energy_array[min_loc + 2], energy_array[min_loc + 1]) + current_time;
            //tmp_double = (pt->time - current_time)*sqrt(energy_array[min_loc + 2] + old_e_right)/sqrt(energy_array[min_loc + 2] + energy_array[min_loc + 1]) + current_time;
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

    int N = 50;
    // if(argc > 1)
    // {
    //     N = strtol(argv[1], NULL,10 );
    // }
    double big_tau = 0.2;//big time step of tau leaping
    const int ratio = int(N/10);//ratio of big step and small step
    double small_tau = big_tau/double(ratio);//small time step

    const int N_thread = 8;
    // double *energy_integration = new double[N_thread];
    //
    // for(int i = 0 ; i < N_thread ; i++) {
    //     energy_integration[i] = 0.0;
    // }

    array<double, N_thread> energy_integration;
    array<thread_info, N_thread> infos;

    for(auto&& item : energy_integration) {
        item = 0;
    }

    int Step = 1000000;

    //omp_lock_t lck;
    //omp_init_lock(&lck);

    #pragma omp parallel num_threads(N_thread)
    {
        int rank = omp_get_thread_num();
        //int size = omp_get_num_threads();

        //cout<<"Hello From "<<rank<<endl;
        // printf("Hello(%d)", rank);
        // printf("World(%d)\n", rank);

        double* energy_array = new double[N+2];
        double *E_avg = new double[N];
        double* last_update = new double[N];
        for(int i = 0; i <N; i++)
        {
            E_avg[i] = 0;
            last_update[i] = 0;
        }

        trng::yarn2 r;
        trng::uniform01_dist<> u;
        r.seed(time(NULL));
        r.split(N_thread, rank);
        energy_array[0] = TL;
        energy_array[N+1] = TR;
        for(int n = 1; n < N+1; n++)
            energy_array[n] = 1;
        interaction* time_array = new interaction[N+1];
        for(int n = 0; n < N+1; n++)
        {
            time_array[n].time = -log(1 - u(r))/rate_function(energy_array[n], energy_array[n+1]);//sqrt(energy_array[n] + energy_array[n+1]);
            time_array[n].location = n;
            time_array[n].left = NULL;
            time_array[n].right = NULL;
        }
        int count = 0;

        interaction** clock_time_in_step = new interaction*[ratio + 1];//each element in the array is the head of a list
        for(int i = 0; i < ratio + 1; i++)
        {
            clock_time_in_step[i] = NULL;
        }

        gettimeofday(&t1,NULL);

        //#pragma omp critical
        //omp_set_lock(&lck);
        // cout<<"OK: "<<rank<<endl;
        // cout<<"OKOK: "<<rank<<endl;
        // cout<<"OKOKOK: "<<rank<<endl;
        // cout<<"OKOKOKOK: "<<rank<<endl;
        for(int out_n = 0; out_n < Step; out_n++)
        {
            big_step_distribute(clock_time_in_step,time_array,N+1,small_tau,ratio,out_n);

            for(int in_n = 0; in_n < ratio; in_n++)
            {

                if(clock_time_in_step[in_n]!= NULL)
                {
                    update(clock_time_in_step, in_n, N, small_tau, ratio, out_n, time_array, energy_array, u, r, count, energy_integration[rank]);
                }
            }

            clock_time_in_step[ratio] = NULL;

        }
        //omp_unset_lock(&lck);

        //cout<<" energy_integration: "<<energy_integration<<endl;
        gettimeofday(&t2, NULL);

        double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u + t2.tv_usec - t1.tv_usec) / 1.e6;

        infos[rank].current_rank = rank;
        infos[rank].count = count;
        infos[rank].running_time = 1000000*delta/double(count);

        delete[] energy_array;
        delete[] E_avg;
        delete[] time_array;
        delete[] clock_time_in_step;
    }

    //omp_destroy_lock(&lck);

    //double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u + t2.tv_usec - t1.tv_usec) / 1.e6;

//    cout << "total CPU time = " << delta <<endl;
    cout<<" N = "<<N <<endl;
    
    auto sum = 0.0;

    for(auto item : energy_integration) {
        cout<<item<<", ";
        sum += item;
    }
    cout<<endl;
    cout<<"Sum: "<<sum<<endl;

    ofstream running_log;
    running_log.open("running_log.txt", ios::trunc);

    for(auto item : infos) {
        running_log<<"rank: "<<item.current_rank<<endl;
        running_log<<"count: "<<item.count<<endl;
        running_log<<"time: "<<item.running_time<<endl;
        running_log<<"================================"<<endl;
    }
    
    for(auto item : energy_integration) {
        running_log<<item<<", ";
    }
    running_log<<endl;
    
    cout<<"Average energy flux = : "<<sum/(N_thread*Step*big_tau)<<endl;

    running_log.close();

    // for(int i = 0 ; i < N_thread ; i++) {
    //     cout<<energy_integration[i]<<", ";
    // }
    //cout<<"seconds per million event is "<< 1000000*delta/double(count)<<endl;
    //myfile<<" N = "<<N <<endl;
    //myfile<< 1000000*delta/double(count)<<"  ";

    // for(auto item : energy_integration) {
    //     cout<<item<<", ";
    // }

    myfile.close();




}
