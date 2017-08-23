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
#include <utility> // std::pair

typedef pair<int, int> Pair;
using namespace std;

double TL = 1.0;
double TR = 2.0;

enum Test 
{
    type1,
    type2,
    type3,
    stop
};

enum Rate_Func 
{
    rf1,
    rf2
};

struct interaction
{
    bool ishorizontal;
    Pair location;
    double time;
    int index;
    int Level;
    interaction* left;
    interaction* right;
};

struct thread_info 
{
    int current_rank = -1;
    double running_time = 0.0;
    int count = 0;
};

Rate_Func rate_func = rf1;

double rate_function(double x, double y) 
{
    if(x >= 0 && y >= 0)
    {
        return (rate_func == rf1) ? sqrt(x*y/(x+y)) : sqrt(x+y);
    }

    else
    {
        cout<<"error, sqrt of negative number! "<<endl;
        return 0;
    }
}

void print_v(double* Array, int size)
{
    for(int i = 0; i < size; i++)
    {
        cout<< Array[i] << "  ";
    }
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
    if( *Link == NULL)
    {
        return;
    }

    if(pt == NULL)
    {
        cout<<"empty pointer !"<<endl;
    }

    else if( pt->left == NULL && pt->right == NULL)
    {
        *Link = NULL;
        return;
    }

    else
    {
        if((*Link) == pt  )
        {
	        if(pt->right == NULL)
            {
                *Link = NULL;
            }
            else
            {
                *Link = pt->right;
            }
        }

        if( pt->right != NULL)
        {
            pt->right->left = pt->left;
        }
        if( pt->left != NULL)
        {
            pt->left->right = pt->right;
        }
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


void big_step_distribute(interaction** &clock_time_in_step, interaction* time_array, const int N, const double small_tau, const int ratio, const int Step, int& error)
//distribute clock times of a big step into vectors that represent small steps. If the clock time is bigger than a big tau, then it is arranged in the right location
{
    for(int i = 0; i < N; i++)
    {
        int tmp = 0;
        if(time_array[i].time > (Step + 1)*ratio*small_tau)
        {
            tmp = ratio;
        }
        else
        {
            tmp = floor((time_array[i].time - ratio*small_tau*Step)/small_tau);
            if(tmp < 0)
            {
                cout<<"Error! tmp = "<<tmp<<" when i = "<<i<<" small tau = "<<small_tau<<endl;
                cout<<"time difference is "<<time_array[i].time - ratio*small_tau*Step<<endl;
                error = 1;
                break;
            }
        }
        time_array[i].Level = tmp;
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
    old_level = pt->Level;
    
    if ( new_time > (Step + 1)*ratio*small_tau )
    {
        new_level = ratio;
    }

    else
    {
        new_level = floor( (new_time - Step*ratio*small_tau)/small_tau );
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
        pt->Level = new_level;
        remove(&(clock_time_in_step[old_level]), pt);
        push_front(&(clock_time_in_step[new_level]), pt);
    }
}



void update(interaction** &clock_time_in_step, const int level, const int N, const double small_tau, const int ratio, const int Step, interaction* time_array, double *energy_array, trng::uniform01_dist<> &u, trng::yarn2 &r, int &count, double &energy_integration, int& error)
//update clock_time_in_step[level]
{
    double next_time = (Step*ratio + level + 1)*small_tau;
    int min_loc = find_min(clock_time_in_step[level]);
    interaction* pt = &time_array[min_loc];
    double current_time = pt->time;
    double previous_energy = 0.0;
    double current_energy = 0.0;
    // cout<<"at level " << level << endl;
    count = 0;

    while(clock_time_in_step[level] != NULL)
    {
        count++;
        
        //Step 1: update min interaction and energy
        double tmp_double;
        double rnd;
        double previous_energy;

        int min_loc_row = (*pt).location.first;
        int min_loc_col = (*pt).location.second;
        
        int x = min_loc_row/2;
        int y = min_loc_col;
    
        if(pt->ishorizontal == 1)
        {
            y = min_loc_col + 1;
        }
            
        interaction* pt_prev = pt;

        do
        {
            rnd = u(r);
        } while(rnd < 1e-15 || rnd > 1 -  1e-15);


        if(pt->ishorizontal == 0)
        {   
            old_e_left = energy_array[x][y];
            old_e_right = energy_array[x][y + 1];

            total_energy = energy_array[x][y] + energy_array[x][y + 1];

            tmp_double = - log(rnd)/rate_function(energy_array[x][y], energy_array[x][y + 1]);
        }
            
        else
        {
            old_e_up = energy_array[x][y];
            old_e_down = energy_array[x + 1][y];

            total_energy = energy_array[x][y] + energy_array[x + 1][y];

            tmp_double = - log(rnd)/rate_function(energy_array[x][y], energy_array[x + 1][y]);
        }

        previous_energy = energy_array[x][y];

        // tmp_double = - log(rnd)/rate_function(energy_array[min_loc], energy_array[min_loc + 1]);

            if(tmp_double > 1e10 || tmp_double < 1e-15)
            {
                cout << "tmp_double = " << tmp_double << endl;
            }
            
        // cout<<"added time = " << tmp_double << endl;
        double tmp_rnd_uni;

        do
        {
            tmp_rnd_uni = u(r);
        }while(tmp_rnd_uni < 1e-15 || tmp_rnd_uni > 1 -  1e-15);
        

        if(pt->ishorizontal == 0)
        {
            if(y == 0)
            {
                total_energy = old_e_right - log(1 - u(r)) * TL;
            }
                
            if(y == M)
            {
                total_energy = old_e_left - log(1 - u(r)) * TR;
            }
                
            if(y != 0)
            {
                energy_array[x][y] = tmp_rnd_uni * total_energy;
            }
            
            if(y != M)
            {
                energy_array[x][y + 1] = (1 - tmp_rnd_uni) * total_energy;
            }
        }
            
        else
        {
            energy_array[x][y] = (1 - tmp_rnd_uni) * total_energy;
            energy_array[x + 1][y] = tmp_rnd_uni * total_energy;
        }

        //previous_energy = energy_array[min_loc];

        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step,  current_time + tmp_double);
        
        current_energy = energy_array[x][y];
        
        energy_integration += (1 - 2*int(min_loc == 0))*(current_energy - previous_energy);
        //cout<<energy_integration<<endl;

        //Step 2: update other interactions
        if(pt->ishorizontal == 0)
        {
            // Left
            if(y > 0)
            {
                pt = &time_array[min_loc - 1];
                tmp_double = (pt->time - current_time)*rate_function(energy_array[x][y - 1], old_e_left)/
                    rate_function(energy_array[x][y - 1], energy_array[x][y]) + current_time;
                move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
                //cout << "vertical - left clock: " << (*pt).index << endl;
            }

            // Right
            if(y < M)
            {
                pt = &time_array[min_loc + 1];
                tmp_double = (pt->time - current_time)*rate_function(energy_array[x][y + 2], old_e_right)/
                    rate_function(energy_array[x][y + 2] + energy_array[x][y + 1]), current_time;
                move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
                //cout << "vertical - right clock: " << (*pt).index << endl;
            }
            
            // Upper
            if(x > 0)
            {
                // Upper-left
                if(y > 0)
                {
                    pt = &time_array[min_loc - M - 1];
                    tmp_double = (pt->time - current_time)*rate_function(energy_array[x - 1][y], old_e_left)/
                        rate_function(energy_array[x - 1][y], energy_array[x][y]) + current_time;
                    move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
                    //cout << "vertical - upper left clock: " << (*pt).index << endl;
                }
                    
                // Upper-right
                if(y < M)
                {
                    pt = &time_array[min_loc - M];
                    tmp_double = (pt->time - current_time)*rate_function(energy_array[x - 1][y + 1], old_e_right)/
                        rate_function(energy_array[x - 1][y + 1], energy_array[x][y + 1]) + current_time;
                    move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
                    //cout << "vertical - upper right clock: " << (*pt).index << endl;
                }
            }
                
            // Lower
            if(x < N - 1)
            {
                // Lower left
                if(y > 0)
                {
                    pt = &time_array[min_loc + M];
                    tmp_double = (pt->time - current_time)*rate_function(energy_array[x + 1][y], old_e_left)/
                        rate_function(energy_array[x + 1][y], energy_array[x][y]) + current_time;
                    move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
                    //cout << "vertical - lower left clock: " << (*pt).index << endl;
                }
                    
                // Lower right
                if(y < M)
                {
                    pt = &time_array[min_loc + M + 1];
                    tmp_double = (pt->time - current_time)*rate_function(energy_array[x + 1][y + 1], old_e_right)/
                        rate_function(energy_array[x + 1][y + 1], energy_array[x][y + 1]) + current_time;
                    move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
                    //cout << "vertical - lower right clock: " << (*pt).index << endl;
                }
            }
        }


        // Case 2: Horizontal
        else
        {           
            // Upper
            if(x > 0)
            {
                pt = &time_array[min_loc - 2 * M - 1];
                tmp_double = (pt->time - current_time)*rate_function(energy_array[x - 1][y], old_e_up)/
                    rate_function(energy_array[x - 1][y] + energy_array[x][y]), current_time;
                move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
                //cout << "horizontal - top clock: " << (*pt).index << endl;
            }
            
            // Lower
            if(x < N - 2)
            {
                pt = &time_array[min_loc + 2 * M + 1];
                tmp_double = (pt->time - current_time)*rate_function(energy_array[x + 2][y], old_e_down)/
                    rate_function(energy_array[x + 2][y] + energy_array[x + 1][y]) + current_time;
                move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
                //cout << "horizontal - bottom clock: " << (*pt).index << endl;
            }
                
            // Left
            if(y >= 0)
            {
                // Left-upper
                if(x >= 0)
                {
                    pt = &time_array[min_loc - M - 1];
                    tmp_double = (pt->time - current_time)*rate_function(energy_array[x][y - 1], old_e_up)/
                        rate_function(energy_array[x][y - 1], energy_array[x][y]) + current_time;     
                    move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
                    //cout << "horizontal - left upper clock: " << (*pt).index << endl;
                }
                    
                // Left-lower
                if(x < N - 1)
                {
                    pt = &time_array[min_loc + M];
                    tmp_double = (pt->time - current_time)*rate_function(energy_array[x + 1][y - 1], old_e_down)/
                        rate_function(energy_array[x + 1][y - 1], energy_array[x + 1][y]) + current_time;
                    move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
                    //cout << "horizontal - left lower clock: " << (*pt).index << endl;
                }
            }
            
            // Right
            if(y < M + 1)
            {
                // Right-upper
                if(x >= 0)
                {
                    pt = &time_array[min_loc - M];
                    tmp_double = (pt->time - current_time)*rate_function(energy_array[x][y + 1], old_e_up)/
                        rate_function(energy_array[x][y + 1], energy_array[x][y]) + current_time;
                    move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
                    //cout << "horizontal - right upper clock: " << (*pt).index << endl;
                }
                    
                // Right-lower
                if(x < N - 1)
                {
                    pt = &time_array[min_loc + M + 1];
                    tmp_double = (pt->time - current_time)*rate_function(energy_array[x + 1][y + 1], old_e_down)/
                        rate_function(energy_array[x + 1][y + 1], energy_array[x + 1][y]) + current_time;
                    move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
                    //cout << "horizontal - right lower clock: " << (*pt).index << endl;
                }
            }
        }

        
        //Step 3: update current time
        if(clock_time_in_step[level] != NULL)
        {
            min_loc = find_min(clock_time_in_step[level]);
            pt = &time_array[min_loc];
            current_time = pt->time;

                if(current_time > next_time + 1e-9 || current_time < next_time - small_tau)
                {
                    cout.precision(17);
                    cout<<"Error ! current_time = "<<current_time<<" next_time  = "<<next_time<<" diff = "<<current_time - next_time<<endl;
                    cout<<"min_loc = "<<min_loc<<" clock_time_in_step[level] = "<<clock_time_in_step[level]->time<<endl;
                    cout<<" level = "<<level<<" count = "<<count<<endl;

                    for(int i = 0; i < ratio + 1; i++)
                    {
                        cout<<"level "<<i<<" ";
                        print_list(clock_time_in_step[i]);
                    }

            		print_v(energy_array, N+2);

            		for(int i = 0; i < N+1; i++)
                    {
                        cout<<time_array[i].time<<" ";
                    }

            		cout<<endl;
            		error = 1;
            		break;
                }
        }
        // else
        // {
        //    current_time = next_time + 1;
        // }
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

    const int N_thread = 8;
    const int Step = 1000000;

    int ratio = -1;
    double small_tau = 0.0;
    const auto number_of_trails = 20;
    auto sum = 0.0;

    array<double, N_thread> energy_integration;
    //array<thread_info, N_thread> infos;

    const int N_value = 11;
    const int N_Normal = 11;
    const int N_Max = 11;

    const int M_value = 11;
    const int M_Normal = 11;
    const int M_Max = 11;

    const double TR_Max = 3.0;
    const double TL_Max = 10.0;

    int N = N_value;
    int M = M_value;

    // if(argc > 1)
    // {
    //     N = strtol(argv[1], NULL,10 );
    // }
    //cout<<"OK";

    data<<"N,TL,TR,current_trail,";

    for(int i = 1 ; i <= N_thread ; i++) 
    {
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
//        case type3: goto test3; break;
//        case stop: goto terminate; break;
    }
    
test1:
    if(iteration == number_of_trails && N >= N_Max) {
        current_test = type2;
        data<<",,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
        data<<",,,,,,,,,,,,,,"<<endl;
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
            data<<",,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
            data<<",,,,,,,,,,,,,,"<<endl;
            cout<<"----------------------"<<endl;
            start_row += (number_of_trails + 2);
            iteration = 0;
            //iteration ++;
        }
        else if(iteration == number_of_trails && N < 10) {
            N ++;
            data<<",,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
            data<<",,,,,,,,,,,,,,"<<endl;
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
        data<<",,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
        data<<",,,,,,,,,,,,,,"<<endl;
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
            data<<",,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
            data<<",,,,,,,,,,,,,,"<<endl;
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
        data<<",,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
        data<<",,,,,,,,,,,,,,"<<endl;
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
            data<<",,,,,,,,,,,,,=STDEV.S("<<col_in_csv<<start_row<<":"<<col_in_csv<<start_row + number_of_trails - 1<<"),"<<endl;
            data<<",,,,,,,,,,,,,,"<<endl;
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
        
        double *E_avg = new double[N];
        double* last_update = new double[N];

        for(int i = 0; i < N; i++)
        {
            E_avg[i] = 0;
            last_update[i] = 0;
        }
            trng::yarn2 r;
	    r.seed(time(NULL));
        trng::uniform01_dist<> u;
        r.split(N_thread, rank);

        double** energy_array = new double*[N];
            

        for(int j = 0; j < N; j++)
        {
            energy_array[j] = new double[M + 2];
        }
            
        for(int i = 0; i < N; i++)
        {
            energy_array[i][0] = TL;
            energy_array[i][M + 1] = TR;
        }
            
        for(int n = 0; n < N; n++)
        {
            for(int m = 1; m < M + 1; m++)
            {
                energy_array[n][m] = 1; 
            }
        }


        interaction* time_array = new interaction[N+1];

        for(int n = 0; n < 2 * N - 1; n++)
        {   
            // If the clock is on the vertical side.
            if(n % 2 == 0)
            {
                for(int m = 0; m < M + 1; m++)
                {
                    time_array[index].location = make_pair(n, m);
                    time_array[index].ishorizontal = 0;
                    time_array[index].time = -log(1 - u(mt))/
                        rate_function(energy_array[n/2][m], energy_array[n/2][m + 1]);
                    time_array[index].left = NULL;
                    time_array[index].right = NULL;
                    time_array[index].index = index;
                    index++;
                }
            }
                
            // If the clock is on the horizontal side.
            else
            {
                for(int m = 0; m < M; m++)
                {
                    time_array[index].location = make_pair(n, m);
                    time_array[index].ishorizontal = 1;
                    time_array[index].time = -log(1 - u(mt))/
                        rate_function(energy_array[(n-1)/2][m], energy_array[(n-1)/2 + 1][m]);
                    time_array[index].left = NULL;
                    time_array[index].right = NULL;
                    time_array[index].index = index;
                    index++;
                }
            }
        }
        
        int count = 0;
        int error = 0;

        interaction** clock_time_in_step = new interaction*[ratio + 1];
        //each element in the array is the head of a list

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
            big_step_distribute(clock_time_in_step, time_array, N + 1 , small_tau, ratio, out_n, error);
            
            if(error == 1)
            {
                cout<<"floor = "<<ratio*small_tau*out_n<<endl;
                for(int i = 0; i < N+1; i++)
                    cout<<time_array[i].time<<"  ";
                cout<<endl;
                break;
            }
            
            
            for(int in_n = 0; in_n < ratio; in_n++)
            {
                
                if(clock_time_in_step[in_n]!= NULL)
                {
		  update(clock_time_in_step, in_n, N, small_tau, ratio, out_n, time_array, energy_array, u, r, count, energy_integration[rank], error);
                }
		if(error == 1)
		  break;
            }
	    if(error == 1)
	      break;
            
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
	delete[] last_update;
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
