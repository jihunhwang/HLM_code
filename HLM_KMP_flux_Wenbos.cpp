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
#include <string>
using namespace std;

double TL = 1.0;
double TR = 2.0;

enum Test {
    type1,
    type2,
    type3,
    stop
};

enum Rate_Func {
    rf1,
    rf2
};

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

Rate_Func rate_func = rf1;

inline double rate_function(double x, double y) {
    
    return (rate_func == rf1) ? sqrt(x*y/(x+y)) : sqrt(x+y);
    
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
    
    if ( new_time >= (Step + 1)*ratio*small_tau )
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
    // cout<<"at level " << level << endl;

    while(current_time < next_time && clock_time_in_step[level] != NU
        LL)
    //while(clock_time_in_step[level] != NULL)
    {
        count++;
        
        //Step 1: update min interaction and energy
        double total_energy = energy_array[min_loc] + energy_array[min_loc + 1];
        double tmp_double;
        double tmp_rate;
        double rnd;

        do{
            rnd = u(r);
        }while(rnd < 1e-16 || rnd > 1 - 1e-16);

        tmp_rate = rate_function(energy_array[min_loc], energy_array[min_loc + 1]);

        if(tmp_rate < 1e-16|| std::isnan(tmp_rate))
        {
	       tmp_rate = 1e-16;
	       cout<<"energy ="<< energy_array[min_loc]<<" "<<energy_array[min_loc + 1]<<endl;
        }

        tmp_double = - log(rnd)/tmp_rate;
        //tmp_double = - log(rnd)/tmp_rate + current_time;

        //        cout<<"added time = " << tmp_double << endl;
        double old_e_left = energy_array[min_loc];
        double old_e_right = energy_array[min_loc + 1];
        double tmp_rnd_uni;

        do{
            tmp_rnd_uni = u(r);
        }while(tmp_rnd_uni < 1e-16 || tmp_rnd_uni > 1 - 1e-16);
        
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

        if(min_loc == 0) 
        {
            current_energy = energy_array[min_loc + 1];
        }
        
        energy_integration += (1 - 2*int(min_loc == 0))*(current_energy - previous_energy);
        //cout<<energy_integration<<endl;

        //Step 2: update other interactions
        if(min_loc > 0)
        {
            pt = &time_array[min_loc - 1];
            tmp_double = (pt->time - current_time)*rate_function(energy_array[min_loc - 1], old_e_left)/
                rate_function(energy_array[min_loc - 1], energy_array[min_loc]) + current_time;

            if(tmp_double > 1e16 || tmp_double < current_time || std::isnan(tmp_double))
            {
                cout<<"tmp_double = "<<tmp_double<<endl;		
                cout<<"energy ="<< energy_array[min_loc - 1]<<" "<<energy_array[min_loc]<<endl;
            }

            if(pt->time < current_time)
            {
                cout<<"Error with time"<<pt->time<<"  "<<current_time<<endl;
            }
            //tmp_double = (pt->time - current_time)*sqrt(energy_array[min_loc - 1] + old_e_left)/sqrt(energy_array[min_loc - 1] + energy_array[min_loc]) + current_time;
            move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
        }

        if(min_loc < N)
        {
            pt = &time_array[min_loc + 1];
            tmp_double = (pt->time - current_time)*rate_function(energy_array[min_loc + 2], old_e_right)/
                rate_function(energy_array[min_loc + 2], energy_array[min_loc + 1]) + current_time;

            if(tmp_double > 1e16 || tmp_double < current_time || std::isnan(tmp_double))
            {
                cout<<"tmp_double = "<<tmp_double<<endl;
        		cout<<"energy ="<< energy_array[min_loc + 1]<<" "<<energy_array[min_loc + 2]<<endl;
        	}

            if(pt->time < current_time)
            {
                cout<<"Error with time"<<pt->time<<"  "<<current_time<<endl;
            }
            //tmp_double = (pt->time - current_time)*sqrt(energy_array[min_loc + 2] + old_e_right)/sqrt(energy_array[min_loc + 2] + energy_array[min_loc + 1]) + current_time;
            move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
            
        }
        
        //Step 3: update current time
        if(clock_time_in_step[level] != NULL)
        {
            min_loc = find_min(clock_time_in_step[level]);
            pt = &time_array[min_loc];
            current_time = pt->time;

            if(current_time > next_time)
            {
                cout << "Error! " << "current_time: " << current_time << " pt->time: " << pt->time << endl;
            }
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
    ofstream data;
    myfile.open("HL_KMP.txt", ios_base::app);
    data.open("data.csv", ios_base::trunc);
    Test current_test = type1;
    int iteration = 0;
    int seed_flag = 0;
    double big_tau = 0.2;//big time step of tau leaping
    const int N_thread = 4;
    const int Step = 100000000;
    int ratio = -1;
    double small_tau = 0.0;
    const auto number_of_trails = 10;
    auto sum = 0.0;
    array<double, N_thread> energy_integration;
    //array<thread_info, N_thread> infos;
    const int N_value = 60;
    const int N_Normal = 10;
    const int N_Max = 60;
    const double TR_Max = 3.0;
    const double TL_Max = 10.0;
    int N = N_value;
    // if(argc > 1)
    // {
    //     N = strtol(argv[1], NULL,10 );
    // }
    //cout<<"OK";
    data<<"N,TL,TR,current_trail,";
    for(int i = 1 ; i <= N_thread ; i++) {
        data<<"thread_"<<i<<",";
    }
    data<<"Mean,std_div"<<endl;
    data<<"Rate Function: sqrt(x*y/(x+y)),,,,,,,,,,,,,,"<<endl;
    data<<"Test_1: Mean energy flux for changing length of the chain. Let N change from 2 to 100.,,,,,,,,,,,,,,"<<endl;
    char col_in_csv = ('A' + 5 + N_thread - 1);
    int start_row = 4;
    //bool printed = false;
    //cout<<col_in_csv<<endl;
second_part:
    if(rate_func == rf2) {
        data<<"Rate Function: sqrt(x+y),,,,,,,,,,,,,,"<<endl;
        data<<"Test_1: Mean energy flux for changing length of the chain. Let N change from 2 to 100.,,,,,,,,,,,,,,"<<endl;
    }
    
check:
    switch (current_test) {
        case type1: goto test1; break;
        case type2: goto terminate; break;
	  //        case type2: goto test2; break;
	  //  case type3: goto test3; break;
	  //  case stop: goto terminate; break;
    }
    
test1:
    if(iteration == number_of_trails && N >= N_Max) {
        current_test = type2;
        data<<",,,,,,,,,,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
	data<<",,,,,,,,,,,,,,,,,,,,,=AVERAGE("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
        data<<",,,,,,,,,,,,,,,,,,,,,,"<<endl;
        data<<"Test_2: Mean energy flux for changing difference of boundary temperatures. Fix TL = 1 and change TR from 1.5 to 10.0.,,,,,,,,,,,,,,"<<endl;
        cout<<"----------------------"<<endl;
        start_row += (number_of_trails + 3);
        iteration = 0;
        N = N_Normal;
        TR = 1.5;
        goto check;
    }
    else {
        if(iteration == number_of_trails && N >= 10) {
            N += 5;
            data<<",,,,,,,,,,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
	      data<<",,,,,,,,,,,,,,,,,,,,,=AVERAGE("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
            data<<",,,,,,,,,,,,,,,,,,,,,,"<<endl;
            cout<<"----------------------"<<endl;
            start_row += (number_of_trails + 2);
            iteration = 0;
            //iteration ++;
        }
        else if(iteration == number_of_trails && N < 10) {
            N ++;
            data<<",,,,,,,,,,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
	    	data<<",,,,,,,,,,,,,,,,,,,,,=AVERAGE("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
            data<<",,,,,,,,,,,,,,,,,,,,,,"<<endl;
            cout<<"----------------------"<<endl;
            start_row += (number_of_trails + 2);
            iteration = 0;
        }
        iteration ++;
        cout<<" N: "<<N<<" ,current iteration: "<<iteration<<endl;
        goto test_begin;
    }
    
test2:
    if(iteration == number_of_trails && TR >= TR_Max) {
        current_test = type3;
        data<<",,,,,,,,,,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
	data<<",,,,,,,,,,,,,,,,,,,,,=AVERAGE("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
        data<<",,,,,,,,,,,,,,,,,,,,,,"<<endl;
        data<<"Test_3: Mean energy flux for changing temperatures. Change TL from 1 to 100 and let TR be 1.10 times TL.,,,,,,,,,,,,,,"<<endl;
        cout<<"----------------------"<<endl;
        start_row += (number_of_trails + 3);
        iteration = 0;
        TL = 1.0;
        TR = TL * 1.10;
        goto check;
    }
    else {
        if(iteration == number_of_trails) {
            TR += 0.5;
            data<<",,,,,,,,,,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
	    data<<",,,,,,,,,,,,,,,,,,,,,=AVERAGE("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
            data<<",,,,,,,,,,,,,,,,,,,,,,"<<endl;
            cout<<"----------------------"<<endl;
            start_row += (number_of_trails + 2);
            iteration = 0;
        }
        iteration ++;
        cout<<" TR: "<<TR<<" ,current iteration: "<<iteration<<endl;
        goto test_begin;
    }
    //        cout<<"aaaaa";
    //        //goto terminate;
test3:
    if(iteration == number_of_trails && TL >= TL_Max) {
        current_test = stop;
        data<<",,,,,,,,,,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
	data<<",,,,,,,,,,,,,,,,,,,,,=AVERAGE("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
        data<<",,,,,,,,,,,,,,,,,,,,,,"<<endl;
        //data<<"EOF,,,,,,,,,,,,,,"<<endl;
        cout<<"----------------------"<<endl;
        start_row += (number_of_trails + 4);
        iteration = 0;
        TL = 1;
        TR = 2;
        goto check;
    }
    else {
        if(iteration == number_of_trails) {
            TL += 5.0;
            TR = 1.10 * TL;
            data<<",,,,,,,,,,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
	    data<<",,,,,,,,,,,,,,,,,,,,,=AVERAGE("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
            data<<",,,,,,,,,,,,,,,,,,,,,,"<<endl;
            cout<<"----------------------"<<endl;
            start_row += (number_of_trails + 2);
            iteration = 0;
        }
        iteration ++;
        cout<<" TL: "<<TL<<", TR: "<<TR<<", current iteration: "<<iteration<<endl;
        goto test_begin;
    }
    //        cout<<"bbbbb";
    //goto terminate;
    
test_begin:
    
    ratio = (N < 10) ? 2 : int(N/10);//ratio of big step and small step
    small_tau = (N < 10) ? 0.1 : big_tau/double(ratio);//small time step
    
    
    // double *energy_integration = new double[N_thread];
    //
    // for(int i = 0 ; i < N_thread ; i++) {
    //     energy_integration[i] = 0.0;
    // }
    
    
    
    for(auto&& item : energy_integration) {
        item = 0;
    }
    
    
    
    //omp_lock_t lck;
    //omp_init_lock(&lck);
    
#pragma omp parallel num_threads(N_thread)
    {
        //omp_set_lock(&lck);
        
        int rank = omp_get_thread_num();
        
        double* energy_array = new double[N+2];
        double *E_avg = new double[N];
        double* last_update = new double[N];
        for(int i = 0; i <N; i++)
        {
            E_avg[i] = 0;
            last_update[i] = 0;
        }
            trng::yarn2 r;
	    r.seed(time(NULL));
        trng::uniform01_dist<> u;
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
        
        //            double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u + t2.tv_usec - t1.tv_usec) / 1.e6;
        //
        //            infos[rank].current_rank = rank;
        //            infos[rank].count = count;
        //            infos[rank].running_time = 1000000*delta/double(count);
        
        delete[] energy_array;
        delete[] E_avg;
        delete[] time_array;
        delete[] clock_time_in_step;
        //omp_unset_lock(&lck);
        
    }
    
    //ßomp_destroy_lock(&lck);
    
    sum = 0;
    
    data<<N<<","<<TL<<","<<TR<<","<<iteration<<",";
    
    for(auto item : energy_integration) {
        data<<item<<",";
        sum += item;
    }
    
    data<<sum/(N_thread*Step*big_tau)<<","<<endl;
    
    //double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u + t2.tv_usec - t1.tv_usec) / 1.e6;
    
    //    cout << "total CPU time = " << delta <<endl;
    //cout<<" N = "<<N <<endl;
    goto check;
    
    //    test_end:
    //
    //    auto sum = 0.0;
    //
    //    for(auto item : energy_integration) {
    //        cout<<item<<", ";
    //        sum += item;
    //    }
    //    cout<<endl;
    //    cout<<"Sum: "<<sum<<endl;
    //
    //    ofstream running_log;
    //    running_log.open("running_log.txt", ios::trunc);
    //
    //    for(auto item : infos) {
    //        running_log<<"rank: "<<item.current_rank<<endl;
    //        running_log<<"count: "<<item.count<<endl;
    //        running_log<<"time: "<<item.running_time<<endl;
    //        running_log<<"================================"<<endl;
    //    }
    //
    //    for(auto item : energy_integration) {
    //        running_log<<item<<", ";
    //    }
    //    running_log<<endl;
    //
    //    cout<<"Average energy flux = : "<<sum/(N_thread*Step*big_tau)<<endl;
    //
    //    running_log.close();
    
    // for(int i = 0 ; i < N_thread ; i++) {
    //     cout<<energy_integration[i]<<", ";
    // }
    //cout<<"seconds per million event is "<< 1000000*delta/double(count)<<endl;
    //myfile<<" N = "<<N <<endl;
    //myfile<< 1000000*delta/double(count)<<"  ";
    
    // for(auto item : energy_integration) {
    //     cout<<item<<", ";
    // }
terminate:
    if(rate_func == rf1) {
        current_test = type1;
        N = N_value;
        TL = 1.0;
        TR = 2.0;
        rate_func = rf2;
        goto second_part;
    }
    myfile.close();
    data.close();
    return 0;
    
    
}
