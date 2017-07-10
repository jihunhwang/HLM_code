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


const double TL = 1.0;
const double TR = 2.0;

struct interaction
{
    double time;
    
    // int location;
    /*
     Use pairs for location in 2D
    */
    pair<int,int> location;
    
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
pair<int,int> find_min(interaction* Link)
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
        cout<<" location: "<< tmp->location <<" time: " << tmp->time << "  " ;
        tmp = tmp->right;
    }
    cout<<endl;
}*/


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


void big_step_distribute(interaction** &clock_time_in_step, interaction* time_array, 
	const int N, const double small_tau, const int ratio, const int Step)
// Distribute clock times of a big step into vectors that represent small steps. 
// If the clock time is bigger than a big tau, then it is arranged in the right location
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



void move_interaction(interaction** &clock_time_in_step, interaction* pt, 
	const double small_tau, const int ratio, const int Step, const double new_time)
// move the interaction pointed by *pt from old bucket to new bucket
// n_move1: move without relinking pointers
// n_move2: move with relinking pointers
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
    
    
    // cout<<"start to move "<< pt->location << " from " << old_level << " to " << new_level<<endl;
    
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

int index_calculator(const int n_row, const int n_col, int a, int b)
{
	if(a % 2 == 0)
	{
		return (int)(a * (n_row + n_col + 1)/2 + b);
	}
	
	else
	{
		return (int)((n_row + n_col + 1) * (a - 1)/2 + n_col + 1 + b);
	}
}

/*
 Updating clocks: N is number of columns and M is number of rows energy array has.
*/
void update(interaction** &clock_time_in_step, const int level, const int N, const int M, const double small_tau, 
	const int ratio, const int Step, interaction* time_array, double **energy_array, 
	uniform_real_distribution<double> &u, mt19937 &mt, int &count)
//update clock_time_in_step[level]
{
	double next_time = (Step*ratio + level + 1)*small_tau;
    
    pair <int,int> min_loc = find_min(clock_time_in_step[level]);
    
    int min_loc_row = min_loc.first;
    int min_loc_col = min_loc.second;
    
    // interaction* pt = &time_array[min_loc_row][min_loc_col];
    
    int time_index = index_calculator(M, N, min_loc_row, min_loc_col);
    
	int time_index1;
	int time_index2;
    
    interaction *pt = &time_array[time_index];
    
    /*
    if(min_loc_row % 2 == 0)
    {
    	int vertical_index = (int)(min_loc_row*(N + M + 1)/2 + min_loc_col);
		interaction* pt = &time_array[vertical_index];
	}
	
	else
	{
		int horizontal_index = (int)((N + M + 1) * (min_loc_row - 1)/2 + N + 1 + b);
		interaction* pt = &time_array[horizontal_index];
	}
	*/
    
    double current_time = pt->time;
	// cout<<"at level " << level << endl;
	
    while(current_time < next_time)
    {
		/*
		 Step 1: update min interaction and energy
		
		 	Two cases: (1) Clock is on the vertical side.
		               (2) Clock is on the horizontal side.
		
		
		 Step 2: update other interactions
		 
		 	Note that we need to update six clocks, not two.
		
		
		 Step 3: update current time
		 	
		 	It shouldn't be a big issue..?
		*/
		
		
		count++;
		double total_energy;
		
		double tmp_double;
		double tmp_double1;
		double tmp_double2;
		
		double old_e_left;
		double old_e_right;
		double old_e_up;
		double old_e_down;
		
		/* Random variable (Uniform) */
		double tmp_rnd_uni;
		
		interaction* pt1;
		interaction* pt2;
		
		
		
		
		// (1) Clock is on the vertical side:
		if(pt->ishorizontal == 0)
		{
			/*
		 	 Step 1: update min interaction and energy
		 	 
			 When the clock is on the vertical edge, its coordinate is in (2*k,q) form (k, q are integers)
			 Which means that (k, q-1) and (k, q) of energy array will have their energy changed.
			 i.e. Energy will be exchanged horizontally.
			*/
			
			total_energy = energy_array[(int)(min_loc_row/2)][(int)min_loc_col - 1] + energy_array[(int)(min_loc_row/2)][min_loc_col];
			tmp_double = -log(1 - u(mt))/sqrt(total_energy);
			
			// old_e_left = energy_array[min_loc_row][];
			// old_e_right = energy_array[min_loc + 1];
			old_e_left = energy_array[(min_loc_row)/2][min_loc_col - 1];
			old_e_right = energy_array[(min_loc_row)/2][min_loc_col];
			
			
			/* 
			 Random variable (Uniform) 
			*/
        	tmp_rnd_uni = u(mt);
        	
        	
        	if(min_loc_col == 0)
        	{
        		total_energy = old_e_right - log(u(mt))*TL;
			}
			
			if(min_loc_col == N)
        	{
            	total_energy = old_e_left - log(u(mt))*TR;
        	}
        	
        	if(min_loc_col != 0)
        	{
            	energy_array[(min_loc_row)/2][min_loc_col - 1] = tmp_rnd_uni * total_energy;
        	}
        
        	if(min_loc_col != N)
        	{
				energy_array[(min_loc_row)/2][min_loc_col] = (1 - tmp_rnd_uni)*total_energy;
        	}
        	
        	move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, current_time + tmp_double);
        	
        	
        	/*
		 	 Step 2: update other interactions
		 	 
		 	 In general, if the clock on (a,b) goes off, the six clocks that are subject to be updated are
		 	 
			 (a-1, b-1), (a-1, b), (a, b-1), (a, b+1), (a+1, b-1), (a+1, b)
			 
			 However, the points with zero or negative coordinate or larger than N+1 or 2M-1 will be excluded.
			*/
			
			// Update the left-hand-side clock
			if(min_loc_col != 0)
        	{
            	time_index = index_calculator(M, N, min_loc_row, min_loc_col - 1);
				pt = &time_array[time_index];
            	
            	tmp_double = (pt->time - current_time)*sqrt(energy_array[min_loc_row/2][min_loc_col - 1] + old_e_left)/
					sqrt(energy_array[min_loc_row/2][min_loc_col - 1] + energy_array[min_loc_row/2][min_loc_col]) + current_time;
					
            	move_interaction(clock_time_in_step, pt, small_tau,ratio, Step, tmp_double);
        	}
        
        	// Update the right-hand-side clock
        	if(min_loc_col != N)
        	{
            	time_index = index_calculator(M, N, min_loc_row, min_loc_col + 1);
				pt = &time_array[time_index];
            	
            	tmp_double = (pt->time - current_time)*sqrt(energy_array[min_loc_row/2][min_loc_col + 2] + old_e_right)/
					sqrt(energy_array[min_loc_row/2][min_loc_col + 2] + energy_array[min_loc_row/2][min_loc_col + 1]) + current_time;
					
            	move_interaction(clock_time_in_step, pt, small_tau,ratio, Step, tmp_double);
        	}
        	
        	// Update two upper clocks
        	if(min_loc_row != 0)
        	{
        		// Upper-left
				if(min_loc_col != 0)
        		{
        			time_index1 = index_calculator(M, N, min_loc_row - 1, min_loc_col);
        			pt1 = &time_array[time_index1];
        			
        			tmp_double1 = (pt1->time - current_time)*sqrt(energy_array[min_loc_row/2 - 1][min_loc_col] + old_e_right)/
						sqrt(energy_array[min_loc_row/2 - 1][min_loc_col] + energy_array[min_loc_row/2 - 2][min_loc_col]) + current_time;
					
					move_interaction(clock_time_in_step, pt1, small_tau,ratio, Step, tmp_double);
				}
				
				// Upper-right
				if(min_loc_col != N)
				{
					time_index2 = index_calculator(M, N, min_loc_row - 1, min_loc_col + 1);
					pt2 = &time_array[time_index2];
					
					tmp_double2 = (pt2->time - current_time)*sqrt(energy_array[min_loc_row/2 - 1][min_loc_col + 1] + old_e_right)/
						sqrt(energy_array[min_loc_row/2 - 1][min_loc_col + 1] + energy_array[min_loc_row/2 - 2][min_loc_col + 1]) + current_time;
					
					move_interaction(clock_time_in_step, pt2, small_tau,ratio, Step, tmp_double);
				}
			}
        
        	// Update two lower clocks
        	if(min_loc_row != M)
        	{
        		// Lower-left
        		if(min_loc_col != 0)
        		{
        			time_index1 = index_calculator(M, N, min_loc_row + 1, min_loc_col);
        			pt1 = &time_array[time_index1];
        			
        			tmp_double1 = (pt1->time - current_time)*sqrt(energy_array[min_loc_row/2 + 1][min_loc_col] + old_e_right)/
						sqrt(energy_array[min_loc_row/2 + 1][min_loc_col] + energy_array[min_loc_row/2 + 2][min_loc_col]) + current_time;
						
					move_interaction(clock_time_in_step, pt1, small_tau,ratio, Step, tmp_double);
				}
				
				// Lower-right
				if(min_loc_col != N)
				{
					time_index2 = index_calculator(M, N, min_loc_row + 1, min_loc_col + 1);
					pt2 = &time_array[time_index2];
					
					tmp_double2 = (pt2->time - current_time)*sqrt(energy_array[min_loc_row/2 + 1][min_loc_col + 1] + old_e_right)/
						sqrt(energy_array[min_loc_row/2 + 1][min_loc_col + 1] + energy_array[min_loc_row/2 + 2][min_loc_col + 1]) + current_time;
						
					move_interaction(clock_time_in_step, pt2, small_tau,ratio, Step, tmp_double);
				}
			}
        	// So, in general, there should be six clocks updated after the reaction (depends on the location).
        	
        	
        	
        	/*
		 	 Step 3: update current time
			*/
        	if(clock_time_in_step[level] != NULL)
        	{
            	min_loc = find_min(clock_time_in_step[level]);
            	
            	min_loc_row = min_loc.first;
            	min_loc_col = min_loc.second;
            	
            	time_index = index_calculator(M, N, min_loc_row, min_loc_col);
            	
            	pt = &time_array[time_index];
            	current_time = pt->time;
        	}
        
        	else
        	{
            	current_time = next_time + 1;
        	}
		}
		
		
		
		// (2) Clock is on the horizontal side:
		else
		{
			/*
		 	 Step 1: update min interaction and energy
		 	 
		 	 
			 When the clock is on the horizontal edge, its coordinate is in (2*k + 1,q) form (k, q are integers)
			 Which means that (k, q) and (k + 1, q) of energy array will have their energy changed.
			 i.e. Energy will be exchanged vertically.
			*/
			
			total_energy = energy_array[(int)((min_loc_row - 1)/2)][(int)min_loc_col]
				 + energy_array[(int)((min_loc_row - 1)/2 + 1)][(int)min_loc_col];
				 
			tmp_double = -log(1 - u(mt))/sqrt(total_energy);
			
			old_e_up = energy_array[(int)((min_loc_row - 1)/2)][(int)min_loc_col];
			old_e_down = energy_array[(int)((min_loc_row - 1)/2 + 1)][(int)min_loc_col];
			
			
			/* 
			 Random variable (Uniform) 
			*/
        	tmp_rnd_uni = u(mt);
        	
        	/*
        	if(min_loc_col == 0)
        	{
        		total_energy = old_e_right - log(u(mt))*TL;
			}
			
			if(min_loc_col == N)
        	{
            	total_energy = old_e_left - log(u(mt))*TR;
        	}
        	
        	if(min_loc_col != 0)
        	{
            	energy_array[min_loc_row][min_loc_col - 1] = tmp_rnd_uni * total_energy;
        	}
        
        	if(min_loc_col != N)
        	{
				energy_array[min_loc_row][min_loc_col] = (1 - tmp_rnd_uni) * total_energy;
        	}
        	*/
        	
        	energy_array[(int)((min_loc_row - 1)/2)][(int)min_loc_col] = tmp_rnd_uni * total_energy;
        	energy_array[(int)((min_loc_row - 1)/2 + 1)][(int)min_loc_col] = (1 - tmp_rnd_uni) * total_energy;
        	
        	move_interaction(clock_time_in_step, pt,small_tau, ratio, Step, current_time + tmp_double);
        	
        	
        	/*
		 	 Step 2: update other interactions
		 	 
		 	 In general, if the clock on (a,b) goes off, the six clocks that are subject to be updated are
		 	 
			 (a-2, b), (a-1, b), (a-1, b+1), (a+1,b), (a+1, b+1), (a+2, b)
			 
			 However, the points with zero or negative coordinate or larger than N+1 or 2M-1 will be excluded.
			*/
			
			// Update the clock located on the top
			if(min_loc_row != 1)
        	{
            	time_index = index_calculator(M, N, min_loc_row - 2, min_loc_col);
				
				pt = &time_array[time_index];
            	
            	tmp_double = (pt->time - current_time)*sqrt(energy_array[(min_loc_row - 1)/2][min_loc_col] + old_e_up)/
					sqrt(energy_array[(min_loc_row - 1)/2][min_loc_col] + energy_array[(min_loc_row - 3)/2][min_loc_col]) + current_time;
					
            	move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
        	}
        
        	// Update the clock located on the bottom
        	if(min_loc_row != M)
        	{
            	time_index = index_calculator(M, N, min_loc_row + 2, min_loc_col);
				
				pt = &time_array[time_index];
            	
            	tmp_double = (pt->time - current_time)*sqrt(energy_array[(min_loc_row + 1)/2][min_loc_col] + old_e_down)/
					sqrt(energy_array[(min_loc_row + 1)/2][min_loc_col] + energy_array[(min_loc_row + 3)/2][min_loc_col]) + current_time;
					
            	move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
        	}
        	
        	// Update two left-hand side clocks
        	if(min_loc_col != 0)
        	{
        		// Left-upper
				if(min_loc_row != 0)
        		{
        			time_index1 = index_calculator(M, N, min_loc_row - 1, min_loc_col);
        			pt1 = &time_array[time_index1];
        			
        			tmp_double1 = (pt1->time - current_time)*sqrt(energy_array[(min_loc_row - 1)/2][min_loc_col] + old_e_up)/
						sqrt(energy_array[(min_loc_row - 1)/2][min_loc_col] + energy_array[(min_loc_row - 1)/2][min_loc_col - 1]) + current_time;
					
					move_interaction(clock_time_in_step, pt1, small_tau,ratio, Step, tmp_double);
				}
				
				// Left-lower
				if(min_loc_row != M)
				{
					time_index2 = index_calculator(M, N, min_loc_row + 1, min_loc_col);
					pt2 = &time_array[time_index2];
					
					tmp_double2 = (pt2->time - current_time)*sqrt(energy_array[(min_loc_row + 1)/2][min_loc_col] + old_e_down)/
						sqrt(energy_array[(min_loc_row + 1)/2][min_loc_col] + energy_array[(min_loc_row + 1)/2][min_loc_col - 1]) + current_time;
					
					move_interaction(clock_time_in_step, pt2, small_tau,ratio, Step, tmp_double);
				}
			}
        
        
        	// Update two right-hand side clocks
        	if(min_loc_col != N)
        	{
        		// right-upper
        		if(min_loc_row != 0)
        		{
        			time_index1 = index_calculator(M, N, min_loc_row - 1, min_loc_col);
        			pt1 = &time_array[time_index1];
        			
        			tmp_double1 = (pt1->time - current_time)*sqrt(energy_array[(min_loc_row - 1)/2][min_loc_col] + old_e_up)/
						sqrt(energy_array[(min_loc_row - 1)/2][min_loc_col] + energy_array[(min_loc_row - 1)/2][min_loc_col + 1]) + current_time;
						
					move_interaction(clock_time_in_step, pt1, small_tau,ratio, Step, tmp_double);
				}
				
				// right-lower
				if(min_loc_row != M)
				{
					time_index2 = index_calculator(M, N, min_loc_row + 1, min_loc_col);
					pt2 = &time_array[time_index2];
					
					tmp_double2 = (pt2->time - current_time)*sqrt(energy_array[(min_loc_row + 1)/2][min_loc_col] + old_e_down)/
						sqrt(energy_array[(min_loc_row + 1)/2][min_loc_col] + energy_array[(min_loc_row + 1)/2][min_loc_col + 1]) + current_time;
						
					move_interaction(clock_time_in_step, pt2, small_tau,ratio, Step, tmp_double);
				}
			}
        	// So, in general, there should be six clocks updated after the reaction.
        	
        	
        	/*
		 	 Step 3: update current time
			*/
        	if(clock_time_in_step[level] != NULL)
        	{
            	min_loc = find_min(clock_time_in_step[level]);
            	
            	min_loc_row = min_loc.first;
            	min_loc_col = min_loc.second;
            	
            	time_index = index_calculator(M, N, min_loc_row, min_loc_col);
            	
            	pt = &time_array[time_index];
            	current_time = pt->time;
        	}
        
        	else
        	{
            	current_time = next_time + 1;
        	}
        
        
    }
    
}
}


/*
 Create the time array (an array of clocks (interactions))
*/
interaction* make_time_array(int N, int M)
{
	/*
	N represents number of rows
	M represents number of columns
	*/
	
	// Total 2NM + N - M clocks
	// num_clocks = 2N * M + N - M;
	
	// First, make a 2D array of the clocks
	interaction twod_array[2*N - 1][M + 1];
	
	
	// Initialize 2D array first
	for(int i = 0; i < 2*N - 1; i++)
	{
		for(int j = 0; j < M + 1; j++)
		{
			// Vertical
			if(i % 2 == 0)
			{
				twod_array[i][j].location = make_pair(i, j);
				twod_array[i][j].ishorizontal = 0;
			}
			
			// Horizontal
			else
			{
				// Since the clocks on horizontal line can't have M index
				if(j == M)
				{
					twod_array[i][j].location = make_pair(-1,-1);
				}
				
				else
				{
					twod_array[i][j].location = make_pair(i, j);
					twod_array[i][j].ishorizontal = 1;
				}
				
				/*
				twod_array[i][j].location = make_pair(i, j);
				twod_array[i][j].ishorizontal = 1;
				*/
			}
		}
	}
	
	// cout << twod_array[0][0] << endl;
	
	// Now, let's map 2D array to 1D array
	
	// interaction return_array[2 * N * M + N - M];
	interaction *return_array;
	
	int index = 0;
	
	for(int i = 0; i < 2*N - 1; i++)
	{
		for(int j = 0; j < M + 1; j++)
		{
			// interaction *element = &twod_array[i][j];
			
			if(twod_array[i][j].location.first == -1)
			{
				// pass;
			}
			
			else
			{
				return_array[index] = twod_array[i][j];
				return_array[index].index = index;
				index++;
			}
			
			/*
			return_array[index] = element;
			return_array[index].index = index;
			index++;*/
		}
	}
	
	return return_array;
}


int main(int argc, char** argv)
{
    struct timeval t1, t2;
    ofstream myfile;
    myfile.open("HL_KMP.txt", ios_base::app);
    
    // N represents the number of rows
    const int N = 100; 
    
    // M represents the number of columns
    const int M = 100;
    
    /*
    if(argc > 1)
    {
        N = strtol(argv[1], NULL,10 );
    }*/
    
    double big_tau = 0.2; // big time step of tau leaping
    const int ratio = int(N/10); // ratio of big step and small step
    
    /*
     Remember. N has to be larger than 10 no matter what
     Otherwise, small_tau goes to infinity
    */
    double small_tau = big_tau/double(ratio); // small time step
    
    
    // double* energy_array = new double[N+2];
    
    /*
    What are those two attributes for?
    */
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
    
    
    // energy_array[0] = TL;
    // energy_array[N+1] = TR;
    
    /*
     Here, we are going to define the Energy array.
     Since there are N rows and M columns, 
	 the Energy array would have N+2 rows and M+2 columns.
	 
	 Indeed, 0th column, N+1th column will each have energy TL and TR respectively,
	 and 0th row and M+1th row will have NULL (since the horizontal boundaries are insulated)
	 
	 Rest of the cells will have energy 1.
    */
    // double energy_array[N+2][M+2];
    double** energy_array;
    
    for(int i = 0; i < N + 2; i++)
    {
    	energy_array[i][0] = TL;
    	energy_array[i][M + 1] = TR;
	}
	
	for(int j = 0; j < M + 2; j++)
	{
		energy_array[0][j] = 0;
		energy_array[N + 1][j] = 0;
	}
	
    
    for(int n = 1; n < N + 1; n++)
    {
    	for(int m = 1; m < M + 1; m++)
    	{
    		energy_array[n][m] = 1;	
		}
        // energy_array[n] = 1;
	}
	
    // interaction* time_array = new interaction[N+1];
    
    /*
     Initialize the time_array (which works as the clock array in this case)
    */
    interaction* time_array = make_time_array(N, M);
    
    cout << time_array << endl;
    
    int num_clocks = 2 * N * M + N - M; // Number of clocks in general
    
    for(int m = 0; m < M + 1; m++)
    {	
    	for(int n = 0; n < N + 1; n++)
    	{
			int time_index = index_calculator(N, M, n, m);
			
			time_array[time_index].time = -log(1 - u(mt))/sqrt(energy_array[n][m] + energy_array[n + 1][m]);
        	time_array[time_index].location.first = n;
        	time_array[time_index].location.second = m;
        	time_array[time_index].left = NULL;
        	time_array[time_index].right = NULL;
    	}
	}
    
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
	// cout<<"at Step "<<out_n<<endl;
	
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
                update(clock_time_in_step, in_n, M, N, small_tau, ratio, out_n, time_array, energy_array, u, mt, count);
            }
            
		// print_v(energy_array,N+2);
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
    
    
    for(int i = 0; i < N; i++)
    {
        E_avg[i]= E_avg[i]/(Step*big_tau);
    }
    
    print_v(E_avg, N);
    
    
    
    gettimeofday(&t2, NULL);
    delete energy_array;
    delete[] E_avg;
    delete[] time_array;
    delete[] clock_time_in_step;

    double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
                    t2.tv_usec - t1.tv_usec) / 1.e6;
    
	// cout << "total CPU time = " << delta <<endl;
	
    cout<<" N = "<<N <<endl;
    cout<<" M = "<<M <<endl;
    
    // cout<<"seconds per million event is "<< 1000000*delta/double(count)<<endl;
    // myfile<<" N = "<<N <<endl;
    
    myfile<< 1000000*delta/double(count)<<"  ";
    myfile.close();

 
 
    
}

