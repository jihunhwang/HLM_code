#include <iostream>
#include <fstream>
#include <random>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <list>

/*
We will be using the pair for the location
*/
#include <stdio.h>
#include <utility> // std::pair

using namespace std;

typedef pair<int, int> Pair;

const double TL = 1.0;
const double TR = 2.0;

struct interaction
{
    double time;
    
    /*
     Use pairs for location in 2D
    */
    Pair location;
    
    // saving the index of the time array just for convenience
    int index;
    
    /*
	 checking if a clock is located on the horizontal or vertical side
	  -> in fact, if the x-coordinate of the clock is divisible by 2, it is vertical.
	*/
    bool ishorizontal;  /* stores 1 (True) if horizontal, 0 (False) if vertical */
    
    /*
	 Two pointers for the doubly linkedlist in each bucket 
	*/
    interaction* left;
    interaction* right;
    
};



void print_v(double* Array, int size)
{
    for(int i = 0; i < size; i++)
        cout<< Array[i] << "  ";
    cout<<endl;
}



/*
 Find a node (interaction) in a bucket that holds the smallest time
*/
int find_min(interaction* Link)
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
    
    return pt->index;
    
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
	{
		return;
	}
	
    
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
        /*
         Check if pt is the first element of the Link
        */
		if((*Link) == pt )
		{
            *Link = pt->right;
		}
            
        if( pt->right != NULL)
        {
            pt->right->left = pt->left;
		}
		
		/*
         Check if pt is the last element of the Link
        */
        if( pt->left != NULL)
        {
            pt->left->right = pt->right;
		}
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
    
	/*
	 Cut pt's left link first
	*/
	pt->left = NULL;
	
	/*
     if the Link is empty, just put pt into the Link
    */
    if(*Link == NULL)
    {
        pt->right = NULL;
        *Link = pt;
    }
    
    /*
     Otherwise, connect the pt with the head of the Link, 
	 and let pt be the first element of the Link
    */
    else
    {
        pt->right = *Link;
        (*Link)->left = pt;
        *Link = pt;
    }
}


/*
void print_list(interaction* Link)
{
    interaction* tmp = Link;
    while(tmp!= NULL)
    {
        cout<<" location: "<< tmp->location.first << "," << tmp->location.second <<" time: " << tmp->time << "  " ;
        tmp = tmp->right;
    }
    cout<<endl;
}*/



void big_step_distribute(interaction** &clock_time_in_step, interaction* time_array, 
	const int N, const double small_tau, const int ratio, const int Step)
// Distribute clock times of a big step into vectors that represent small steps. 
// If the clock time is bigger than a big tau, then it is arranged in the right location
{
    //cout << "big_step_distribue begins" << endl;
	
	for(int i = 0; i < N; i++)
    {
        // cout << "00" << endl;
		
		int tmp;
		
        if(time_array[i].time > (Step + 1)*ratio*small_tau)
        {
            //cout << "00" << endl;
			tmp = ratio;
        }
        else
        {
            tmp = int((time_array[i].time - ratio*small_tau*Step)/small_tau);
        }
        //cout << "tmp: " << tmp << endl;
        push_front(&clock_time_in_step[tmp], &time_array[i]);
    }
}



void move_interaction(interaction** &clock_time_in_step, interaction* pt, 
	const double small_tau, const int ratio, const int Step, const double new_time)
// move the interaction pointed by *pt from old bucket to new bucket
// n_move1: move without relinking pointers
// n_move2: move with relinking pointers
{
    //cout << "1" << endl;
	
	double old_time = pt->time;
    int old_level, new_level;
    pt->time = new_time;
    
    //cout << new_time << endl;
    if (old_time > (Step + 1)*ratio*small_tau )
    {
        old_level = ratio;
    }
    //cout << "1" << endl;
    else
    {
        old_level = int( (old_time - Step*ratio*small_tau)/small_tau );
    }
    //cout << "1" << endl;
    
    if ( new_time > (Step + 1)*ratio*small_tau )
    {
        new_level = ratio;
        //cout << "1" << endl;
    }
    
    else
    {
        new_level = int( (new_time - Step*ratio*small_tau)/small_tau );
    }
    //cout << new_level << endl;
    //cout << new_time << endl;
    //cout << Step*ratio*small_tau << endl;
    
    if(old_level == new_level )
    {
        pt->time = new_time;
        //cout << "1" << endl;
    }
    
    else
    {
	    //cout << "5" << endl;
		remove(&(clock_time_in_step[old_level]), pt);
	    //cout << "1" << endl;
	    //cout << new_level << endl;
        push_front(&(clock_time_in_step[new_level]), pt);
        //cout << "1" << endl;
    }
}


/*
 This function makes finding the index of certain clock (a,b) on the time_array easier and faster.
*//*
int index_calculator(const int n_row, const int n_col, int a, int b)
{
	// Is the clock on the vertical side?
	if(a % 2 == 0)
	{
		return (int)(a * (n_row + n_col + 1)/2 + b);
	}
	
	// Is the clock on the horizontal side?
	else
	{
		return (int)((n_row + n_col + 1) * (a - 1)/2 + n_col + 1 + b);
	}
}*/


/*
 Updating clocks: N is number of rows and M is number of columns energy array has.
*/
void update(interaction** &clock_time_in_step, const int level, const int N, const int M, const double small_tau, 
	const int ratio, const int Step, interaction* time_array, double **energy_array, 
	uniform_real_distribution<double> &u, mt19937 &mt, int *count)
//update clock_time_in_step[level]
{
//	cout << "" << endl;
	//cout << "UPDATE FUNCTION STARTS!" << endl;
	
	double next_time = (Step*ratio + level + 1)*small_tau;
    
    int min_loc = find_min(clock_time_in_step[level]);
    
    interaction *pt = &time_array[min_loc];
    
    //cout << "min_loc: " << min_loc << endl;
    
    double current_time = pt->time;
    
	
	int min_loc_row = (*pt).location.first;
	int min_loc_col = (*pt).location.second;
		
	int x = min_loc_row/2;
	int y = min_loc_col;
	
	//interaction *pt1;
	
    while(current_time < next_time)
    {
		//cout << "WHILE LOOP!" << endl;
		
		//cout << "coordinate: " << (pt->location).first << "," << (pt->location).second << endl;
		//cout << "min_loc: " << min_loc << endl;
		//cout << "pt index: " << pt->index << endl;
		
    	//cout << "next_time: " << next_time << endl;
		//cout << "current_time: " << current_time << endl;
		
		(*count)++;
		

		min_loc_row = (*pt).location.first;
		min_loc_col = (*pt).location.second;
		
		x = min_loc_row/2;
		y = min_loc_col;
	
		/*if(pt->ishorizontal == 1)
		{
			y = min_loc_col + 1;
		}*/
		
		double total_energy, tmp_double;
		double old_e_left, old_e_right, old_e_up, old_e_down;
		
		double tmp_rnd_uni = u(mt);
		
		//cout << "1" << endl;
		
		if(pt->ishorizontal == 0)
		{	
        	//cout << "1" << endl;
			old_e_left = energy_array[x][y];
			old_e_right = energy_array[x][y + 1];
			//cout << "1" << endl;
		}
		
		else
		{
			//cout << "2" << endl;
			old_e_up = energy_array[x][y];
			old_e_down = energy_array[x + 1][y];
		}
        
        //cout << "1" << endl;
		if(pt->ishorizontal == 0)
		{
			//cout << "1" << endl;
			total_energy = energy_array[x][y] + energy_array[x][y + 1];
			//cout << "1" << endl;
		}
		
		else
		{
			total_energy = energy_array[x][y] + energy_array[x + 1][y];
		}
		
		
		tmp_double = -log(1 - u(mt))/sqrt(total_energy);
        //cout << "1" << endl;	

        if(pt->ishorizontal == 0)
        {
        	if(y == 0)
        	{
        		total_energy = old_e_right - log(1 - u(mt)) * TL;
			}
			
			if(y == M)
        	{
            	total_energy = old_e_left - log(1 - u(mt)) * TR;
        	}
        	
        	if(y != 0)
        	{
            	energy_array[x][y] = tmp_rnd_uni * total_energy;
        	}
        
        	if(y != M)
        	{
				energy_array[x][y + 1] = (1 - tmp_rnd_uni) * total_energy;
        	}
        	//cout << "1" << endl;
		}
        
        else
		{
			energy_array[x][y] = (1 - tmp_rnd_uni) * total_energy;
            energy_array[x + 1][y] = tmp_rnd_uni * total_energy;
		}	

        move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, current_time + tmp_double);
		//Part 2;
        
        // Vertical
        if(pt->ishorizontal == 0)
		{
			// Left
			if(y > 0)
        	{
				pt = &time_array[min_loc - 1];
				tmp_double = (pt->time - current_time)*sqrt(energy_array[x][y - 1] + old_e_left)/sqrt(energy_array[x][y - 1] + energy_array[x][y]) + current_time;
				move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
            	//cout << "vertical - left clock: " << (*pt).index << endl;
        	}

			// Right
        	if(y < M)
        	{
				pt = &time_array[min_loc + 1];
            	tmp_double = (pt->time - current_time)*sqrt(energy_array[x][y + 2] + old_e_right)/sqrt(energy_array[x][y + 2] + energy_array[x][y + 1]) + current_time;
				move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
			//	cout << "vertical - right clock: " << (*pt).index << endl;
        	}
        
        	// Upper
        	if(x > 0)
        	{
        		// Upper-left
				if(y > 0)
        		{
        			pt = &time_array[min_loc - M - 1];
					tmp_double = (pt->time - current_time)*sqrt(energy_array[x - 1][y] + old_e_left)/sqrt(energy_array[x - 1][y] + energy_array[x][y]) + current_time;
					move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
					//cout << "vertical - upper left clock: " << (*pt).index << endl;
				}
				
				// Upper-right
				if(y < M)
				{
					pt = &time_array[min_loc - M];
					tmp_double = (pt->time - current_time)*sqrt(energy_array[x - 1][y + 1] + old_e_right)/sqrt(energy_array[x - 1][y + 1] + energy_array[x][y + 1]) + current_time;
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
        			tmp_double = (pt->time - current_time)*sqrt(energy_array[x + 1][y] + old_e_left)/sqrt(energy_array[x + 1][y] + energy_array[x][y]) + current_time;
					move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
					//cout << "vertical - lower left clock: " << (*pt).index << endl;
				}
				
				// Lower right
				if(y < M)
				{
					pt = &time_array[min_loc + M + 1];
					tmp_double = (pt->time - current_time)*sqrt(energy_array[x + 1][y + 1] + old_e_right)/sqrt(energy_array[x + 1][y + 1] + energy_array[x][y + 1]) + current_time;
					move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
					//cout << "vertical - lower right clock: " << (*pt).index << endl;
				}
			}
        }

		// horizontal
		else
		{			
			// Upper
			if(x > 0)
        	{
				pt = &time_array[min_loc - 2 * M - 1];
            	tmp_double = (pt->time - current_time)*sqrt(energy_array[x - 1][y] + old_e_up)/sqrt(energy_array[x - 1][y] + energy_array[x][y]) + current_time;
				move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
            	//cout << "horizontal - top clock: " << (*pt).index << endl;
        	}
        

        	// Lower
        	if(x < N - 2)
        	{
				pt = &time_array[min_loc + 2 * M + 1];
            	tmp_double = (pt->time - current_time)*sqrt(energy_array[x + 2][y] + old_e_down)/sqrt(energy_array[x + 2][y] + energy_array[x + 1][y]) + current_time;
				move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
            	//cout << "horizontal - bottom clock: " << (*pt).index << endl;
        	}
        	
        	
        	// Left
        	if(y > 0)
        	{
        		// Left-upper
        		if(x >= 0)
        		{
        			pt = &time_array[min_loc - M - 1];
        			tmp_double = (pt->time - current_time)*sqrt(energy_array[x][y - 1] + old_e_up)/sqrt(energy_array[x][y - 1]+ energy_array[x][y]) + current_time;		
					move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
					//cout << "horizontal - left upper clock: " << (*pt).index << endl;
				}
				
				// Left-lower
				if(x < N - 1)
        		{
					pt = &time_array[min_loc + M];
					tmp_double = (pt->time - current_time)*sqrt(energy_array[x + 1][y - 1] + old_e_down)/sqrt(energy_array[x + 1][y - 1]+ energy_array[x + 1][y]) + current_time;
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
        			tmp_double = (pt->time - current_time)*sqrt(energy_array[x][y + 1] + old_e_up)/sqrt(energy_array[x][y + 1] + energy_array[x][y]) + current_time;
					move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
					//cout << "horizontal - right upper clock: " << (*pt).index << endl;
				}
				
				// Right-lower
				if(x < N - 1)
				{
					pt = &time_array[min_loc + M + 1];
					tmp_double = (pt->time - current_time)*sqrt(energy_array[x + 1][y + 1] + old_e_down)/sqrt(energy_array[x + 1][y + 1] + energy_array[x + 1][y]) + current_time;
					move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
					//cout << "horizontal - right lower clock: " << (*pt).index << endl;
				}
			}
		}

			
        /*
		 Step 3: update current time
		*/
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


int main(int argc, char *argv[])
{
    struct timeval t1, t2;
    ofstream myfile;
    myfile.open("HL_KMP.txt", ios_base::app);
    
    // N represents the number of rows
    int N = 100; 
    
    // M represents the number of columns
    int M = 100;
    
    /*
    if(argc > 1)
    {
        N = strtol(argv[1], NULL,10 );
        M = strtol(argv[2], NULL,10 );
    }*/
    
    int num_clocks = 2 * N * M + N - M;
    
    double big_tau = 0.2; // big time step of tau leaping
    const int ratio = int(num_clocks/10); // ratio of big step and small step
    
    /*
     Remember. N has to be larger than 10 no matter what
     Otherwise, small_tau goes to infinity
    */
    double small_tau = big_tau/double(ratio); // small time step
    
    
    double *E_avg = new double[N * M];
    double* last_update = new double[M * N];
    
    
    for(int i = 0; i < N * M; i++)
    {
        E_avg[i] = 0;
        last_update[i] = 0;
    }
    
    random_device rd;
    
    mt19937 mt(rd());
    uniform_real_distribution<double> u(0,1);
    
    
    // energy_array[0] = TL;
    // energy_array[N+1] = TR;
    
    
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
        // energy_array[n] = 1;
	}
	
    
    // Initialize the time_array (which works as the clock array in this case)
    // int num_clocks = 2 * N * M + N - M; // Number of clocks in general
    //interaction time_array[num_clocks];
    interaction* time_array = new interaction[num_clocks];

	int index = 0;
	
    for(int n = 0; n <= 2 * N - 2; n++)
    {	
    	// time_array is vertical
    	if(n % 2 == 0)
    	{
    		for(int m = 0; m < M + 1; m++)
    		{
            	time_array[index].location = make_pair(n, m);
				time_array[index].ishorizontal = 0;
				time_array[index].time = -log(1 - u(mt))/sqrt(energy_array[n/2][m] + energy_array[n/2][m + 1]);
        		time_array[index].left = NULL;
        		time_array[index].right = NULL;
        		time_array[index].index = index;
               	index++;
			}
		}
		
        // time_array is horizontal
        else
        {
            for(int m = 0; m < M; m++)
            {
            	time_array[index].location = make_pair(n, m);
				time_array[index].ishorizontal = 1;
				time_array[index].time = -log(1 - u(mt))/sqrt(energy_array[(n-1)/2][m] + energy_array[(n-1)/2 + 1][m]);
        		time_array[index].left = NULL;
        		time_array[index].right = NULL;
        		time_array[index].index = index;
                index++;
			}
			
        }
    }
    
    // cout << "index: " << index << endl;
    
    int count = 0;
    
	// each element in the array is the head of a list
    interaction** clock_time_in_step = new interaction*[ratio + 1];
    
    for(int i = 0; i < ratio + 1; i++)
    {
        clock_time_in_step[i] = NULL;
    }
    
    gettimeofday(&t1,NULL);
    
	// big_step_distribute(clock_time_in_step,time_array,N+1,small_tau,ratio,0);

    
    int Step = 50;
	   
    for(int out_n = 0; out_n < Step; out_n++)
    {
	 	//cout << "" << endl;
		//cout<< "At Step "<< out_n << endl;
	
        big_step_distribute(clock_time_in_step, time_array, num_clocks, small_tau, ratio, out_n);
        
        for(int in_n = 0; in_n < ratio; in_n++)
        {
			if(clock_time_in_step[in_n]!= NULL)
            {
				//cout<< "" << endl;
				//cout<< "in_n: " << in_n << endl;
				update(clock_time_in_step, in_n, N, M, small_tau, ratio, out_n, time_array, energy_array, u, mt, &count);
            }
		// print_v(energy_array,N+2);
        }
        clock_time_in_step[ratio] = NULL;
    }
    
    gettimeofday(&t2, NULL);

	//cout << "" << endl;
	cout<<"Done! " <<endl;
    
    delete[] energy_array;
    delete[] E_avg;
    delete[] time_array;
    delete[] clock_time_in_step;
    
    double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
                    t2.tv_usec - t1.tv_usec) / 1.e6;
	// cout << "total CPU time = " << delta <<endl;
    cout<<" N = "<< N <<endl;
    cout<<" M = "<< M <<endl;
    
    // cout<<"seconds per million event is "<< 1000000*delta/double(count)<<endl;
    // myfile<<" N = "<<N <<endl;
    
    myfile<< 1000000*delta/double(count)<<"  ";
    myfile.close();

}
