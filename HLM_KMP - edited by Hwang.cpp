	#include <iostream>
	#include <fstream>
	#include <random>
	#include <stdlib.h>
	#include <math.h>
	#include <sys/time.h>
	#include <time.h>
	#include <list>

	// We will be using the pair for the location
	#include <stdio.h>
	#include <utility> // std::pair

	using namespace std;
	using std::vector;

	typedef pair<int, int> Pair;
	// typedef pair<double, double> d_Pair;
	typedef list<double> List;

	vector<int> arr;

	const double TL = 1.0;
	const double TR = 2.0;

	struct interaction
	{
	    double time;
	    
	    // Use pairs for location in 2D
	    Pair location;
	    
	    
		// Saving the index of the time array just for convenience
	    int index;
	    
		// Checking if a clock is located on the horizontal or vertical side
		// -> in fact, if the x-coordinate of the clock is divisible by 2, it is vertical.
	    bool ishorizontal;  // stores 1 (True) if horizontal, 0 (False) if vertical
	    
		// Two pointers for the doubly linkedlist in each bucket
	    interaction* left;
	    interaction* right;
	};


	struct profile
	{
		double prev_time;
		double integral_energy;
		//d_Pair* prevtime_energy;
		Pair location;
	};



	void print_v(double* Array, int size)
	{
	    for(int i = 0; i < size; i++)
	        cout<< Array[i] << "  ";
	    cout<<endl;
	}


	// Find a node (interaction) in a bucket that holds the smallest time
	int find_min(interaction* Link)
	{
		double tmp = Link->time;
		
	    // Initializing the starting point
	    interaction* pt = Link;
	    interaction* tmp_pt = Link;
	   
	    //Keep moving right until it reaches TR
	    while(tmp_pt != NULL)
	    {
	        //Check if tmp is smaller than tmp_pt
			if(tmp_pt->time < tmp)
	        {
				tmp = tmp_pt->time;
	        	pt = tmp_pt;
	        }

	        tmp_pt = tmp_pt->right;
	    }
	    return pt->index;
	}


	//remove pt from link but keep pt for future use
	void remove(interaction** Link, interaction* pt)
	{
	    // Note that interaction** represents the buckets
		

		// Check if the array in a bucket is empty or the point is empty
		if( *Link == NULL || pt == NULL)
		{
			return;
		}
		
	    // Check if the bucket has only one interaction
		else if( pt->left == NULL && pt->right == NULL)
	    {
	        *Link = NULL;
	        return;
	    }
	    
	    else
	    {
	        // Check if pt is the first element of the Link
			if((*Link) == pt )
			{
	            *Link = pt->right;
			}
	            
	        if( pt->right != NULL)
	        {
	            pt->right->left = pt->left;
			}
			
	        // Check if pt is the last element of the Link
	        if( pt->left != NULL)
	        {
	            pt->left->right = pt->right;
			}
	    }
	    
	    // Now, isolate the pt that we just extracted
	    pt->left = NULL;
	    pt->right = NULL;
	}



	void push_front(interaction** Link, interaction* pt)
	//push the interaction pointed by pt into the front of the lise
	{
		// Cut pt's left link first
		pt->left = NULL;
		
	    // if the Link is empty, just put pt into the Link
	    if(*Link == NULL)
	    {
	        pt->right = NULL;
	        *Link = pt;
	    }

	    // Otherwise, connect the pt with the head of the Link, and let pt be the first element of the Link
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
	        //cout<<" location: "<< tmp->location.first << "," << tmp->location.second <<" time: " << tmp->time << "  " ;
	        tmp = tmp->right;
	    }

	    cout<<endl;
	}



	void big_step_distribute(interaction** &clock_time_in_step, interaction* time_array, 
		const int N, const double small_tau, const int ratio, const int Step)
	// Distribute clock times of a big step into vectors that represent small steps. 
	// If the clock time is bigger than a big tau, then it is arranged in the right location
	{
	    //cout << "big_step_distribue begins" << endl;
		
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
		double old_time = pt->time;
	    int old_level, new_level;
	    pt->time = new_time;
	    
	    //cout << new_time << endl;
	    
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



	// This function makes finding the index of certain clock (a,b) on the time_array easier and faster.
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
	}


	// Updating clocks: N is number of rows and M is number of columns energy array has.
	void update(interaction** &clock_time_in_step, const int level, const int N, const int M, const double small_tau, 
		const int ratio, const int Step, interaction* time_array, double **energy_array, 
		uniform_real_distribution<double> &u, mt19937 &mt, int *count, profile** E_avg_profile)
	//update clock_time_in_step[level]
	{
		cout << "" << endl;
		cout << "UPDATE FUNCTION STARTS!" << endl;
		
		double next_time = (Step*ratio + level + 1)*small_tau;
	    int min_loc = find_min(clock_time_in_step[level]);
	    interaction *pt = &time_array[min_loc];
	    
	    cout << "min_loc: " << min_loc << endl;
	    
	    double current_time = pt->time;
	    
		int min_loc_row = (*pt).location.first;
		int min_loc_col = (*pt).location.second;
			
		int x = min_loc_row/2;
		int y = min_loc_col;

	    double previous_time1, previous_time2, previous_energy1, previous_energy2;
	    double time_diff1, time_diff2, integrated1, integrated2, previous_integral1, previous_integral2;
	    double total_energy, tmp_double, curr_t;
	    double old_e_left, old_e_right, old_e_up, old_e_down;
		
		interaction* pt_prev = pt;
		
	    while(current_time < next_time)
	    {
			cout << "WHILE LOOP!" << endl;
			
			cout << "coordinate: " << (pt->location).first << "," << (pt->location).second << endl;
			cout << "min_loc: " << min_loc << endl;
			cout << "pt index: " << pt->index << endl;
			
	    	cout << "next_time: " << next_time << endl;
			cout << "current_time: " << current_time << endl;
			
			(*count)++;
			
			min_loc_row = (*pt).location.first;
			min_loc_col = (*pt).location.second;
			
			x = min_loc_row/2;
			y = min_loc_col;
		
			if(pt->ishorizontal == 1)
			{
				y = min_loc_col + 1;
			}
			
			interaction* pt_prev = pt;
			
			double tmp_rnd_uni = u(mt);
			
			// Step 1: update min interaction and energy
			if(pt->ishorizontal == 0)
			{	
				old_e_left = energy_array[x][y];
				old_e_right = energy_array[x][y + 1];
				
				//if(E_avg_profile[x][y].prev_time == 0)
				//{
					E_avg_profile[x][y].prev_time = (*pt).time;
				//}
				
				//if(E_avg_profile[x][y + 1].prev_time == 0)
				//{
					E_avg_profile[x][y + 1].prev_time = (*pt).time; 
				//}
			}
			
			else
			{
				old_e_up = energy_array[x][y];
				old_e_down = energy_array[x + 1][y];
							
				//if(E_avg_profile[x][y].prev_time == 0)
				//{
					E_avg_profile[x][y].prev_time = (*pt).time;
				//}
				
				//if(E_avg_profile[x][y + 1].prev_time == 0)
				//{
					E_avg_profile[x + 1][y].prev_time = (*pt).time;
				//}
			}
	        

			if(pt->ishorizontal == 0)
			{
				total_energy = energy_array[x][y] + energy_array[x][y + 1];
			}
			
			else
			{
				total_energy = energy_array[x][y] + energy_array[x + 1][y];
			}
			

			tmp_double = -log(1 - u(mt))/sqrt(total_energy);

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
			}
	        
	        else
			{
				energy_array[x][y] = (1 - tmp_rnd_uni) * total_energy;
	            energy_array[x + 1][y] = tmp_rnd_uni * total_energy;
			}


	        move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, current_time + tmp_double);
			

			// Step 2: update other interactions
	        
	        // Case 1: Vertical
	        if(pt->ishorizontal == 0)
			{
				// Left
				if(y > 0)
	        	{
					pt = &time_array[min_loc - 1];
					tmp_double = (pt->time - current_time)*sqrt(energy_array[x][y - 1] + old_e_left)/sqrt(energy_array[x][y - 1] + energy_array[x][y]) + current_time;
					move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
	            	cout << "vertical - left clock: " << (*pt).index << endl;
	        	}

				// Right
	        	if(y < M)
	        	{
					pt = &time_array[min_loc + 1];
	            	tmp_double = (pt->time - current_time)*sqrt(energy_array[x][y + 2] + old_e_right)/sqrt(energy_array[x][y + 2] + energy_array[x][y + 1]) + current_time;
					move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
					cout << "vertical - right clock: " << (*pt).index << endl;
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
						cout << "vertical - upper left clock: " << (*pt).index << endl;
					}
					
					// Upper-right
					if(y < M)
					{
						pt = &time_array[min_loc - M];
						tmp_double = (pt->time - current_time)*sqrt(energy_array[x - 1][y + 1] + old_e_right)/sqrt(energy_array[x - 1][y + 1] + energy_array[x][y + 1]) + current_time;
						move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
						cout << "vertical - upper right clock: " << (*pt).index << endl;
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
						cout << "vertical - lower left clock: " << (*pt).index << endl;
					}
					
					// Lower right
					if(y < M)
					{
						pt = &time_array[min_loc + M + 1];
						tmp_double = (pt->time - current_time)*sqrt(energy_array[x + 1][y + 1] + old_e_right)/sqrt(energy_array[x + 1][y + 1] + energy_array[x][y + 1]) + current_time;
						move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
						cout << "vertical - lower right clock: " << (*pt).index << endl;
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
	            	tmp_double = (pt->time - current_time)*sqrt(energy_array[x - 1][y] + old_e_up)/sqrt(energy_array[x - 1][y] + energy_array[x][y]) + current_time;
					move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
	            	cout << "horizontal - top clock: " << (*pt).index << endl;
	        	}
	        
	        	// Lower
	        	if(x < N - 2)
	        	{
					pt = &time_array[min_loc + 2 * M + 1];
	            	tmp_double = (pt->time - current_time)*sqrt(energy_array[x + 2][y] + old_e_down)/sqrt(energy_array[x + 2][y] + energy_array[x + 1][y]) + current_time;
					move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
	            	cout << "horizontal - bottom clock: " << (*pt).index << endl;
	        	}
	        	
	        	// Left
	        	if(y >= 0)
	        	{
	        		// Left-upper
	        		if(x >= 0)
	        		{
	        			pt = &time_array[min_loc - M - 1];
	        			tmp_double = (pt->time - current_time)*sqrt(energy_array[x][y - 1] + old_e_up)/sqrt(energy_array[x][y - 1]+ energy_array[x][y]) + current_time;		
						move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
						cout << "horizontal - left upper clock: " << (*pt).index << endl;
					}
					
					// Left-lower
					if(x < N - 1)
	        		{
						pt = &time_array[min_loc + M];
						tmp_double = (pt->time - current_time)*sqrt(energy_array[x + 1][y - 1] + old_e_down)/sqrt(energy_array[x + 1][y - 1]+ energy_array[x + 1][y]) + current_time;
						move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
						cout << "horizontal - left lower clock: " << (*pt).index << endl;
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
						cout << "horizontal - right upper clock: " << (*pt).index << endl;
					}
					
					// Right-lower
					if(x < N - 1)
					{
						pt = &time_array[min_loc + M + 1];
						tmp_double = (pt->time - current_time)*sqrt(energy_array[x + 1][y + 1] + old_e_down)/sqrt(energy_array[x + 1][y + 1] + energy_array[x + 1][y]) + current_time;
						move_interaction(clock_time_in_step, pt, small_tau, ratio, Step, tmp_double);
						cout << "horizontal - right lower clock: " << (*pt).index << endl;
					}
				}
			}

			// interaction* pt_prev = pt;
				
			// Step 3: update current time
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
		


	        // Step 4: update the energy average profile
	        if(pt_prev->ishorizontal == 0)
	        {
	        	// First, calculate the difference of previous time and current time
	            time_diff1 = pt->time - E_avg_profile[x][y].prev_time;
	            time_diff2 = pt->time - E_avg_profile[x][y + 1].prev_time;
	            //time_diff1 = current_time - E_avg_profile[x][y].prev_time;
	            //time_diff2 = current_time - E_avg_profile[x][y + 1].prev_time;
	            
	            // Then find the previous energies that has to be integrated over time
	            previous_energy1 = old_e_left;
	            previous_energy2 = old_e_right;
	            
	            // Next, calculate the integral
	            integrated1 = time_diff1 * previous_energy1;
	            integrated2 = time_diff2 * previous_energy2;
	            
	            // Last, add the integrals with the previous integrals
	            E_avg_profile[x][y].integral_energy += integrated1;
	            E_avg_profile[x][y + 1].integral_energy += integrated2;
	            
	            // Don't forget to update prev_time
	        	//E_avg_profile[x][y].prev_time = current_time;
	            //E_avg_profile[x][y + 1].prev_time = current_time;
	            //E_avg_profile[x][y].prev_time = pt->time;
	            //E_avg_profile[x][y + 1].prev_time = pt->time;
	            
	            // Done!
	        }

	        else
	        {
	        	// First, calculate the difference of previous time and current time
	            time_diff1 = pt->time - E_avg_profile[x][y].prev_time;
	            time_diff2 = pt->time - E_avg_profile[x + 1][y].prev_time;
	            //time_diff1 = current_time - E_avg_profile[x][y].prev_time;
	            //time_diff2 = current_time - E_avg_profile[x + 1][y].prev_time;
	            
	            // Then find the previous energies that has to be integrated over time
	            previous_energy1 = old_e_up;
	            previous_energy2 = old_e_down;
	            
	            // Next, calculate the integral
	            integrated1 = time_diff1 * previous_energy1;
	            integrated2 = time_diff2 * previous_energy2;
	            
	            // Last, add the integrals with the previous integrals
	            E_avg_profile[x][y].integral_energy += integrated1;
	            E_avg_profile[x + 1][y].integral_energy += integrated2;
	            
	            // Don't forget to update prev_time
	            //E_avg_profile[x][y].prev_time = current_time;
	            //E_avg_profile[x + 1][y].prev_time = current_time;
	            //E_avg_profile[x][y].prev_time = pt->time;
	            //E_avg_profile[x + 1][y].prev_time = pt->time;
	            
	            // Done!
	        }
		
		}
	}


	int main(int argc, char *argv[])
	{
	    struct timeval t1, t2;
	    ofstream myfile;
	    ofstream myprofile;
	    myfile.open("HL_KMP.txt", ios_base::app);
	    myprofile.open("KMP_Profile.txt", ios_base::app);
	    
	    // N represents the number of rows
	    int N = 11;
	    
	    // M represents the number of columns
	    int M = 11;
	    
	    
	    if(argc > 1)
	    {
	        N = strtol(argv[1], NULL,10 );
	        M = strtol(argv[2], NULL,10 );
	    }
	    
	    // Number of clocks (interactions) in N by M energy array.
	    int num_clocks = 2 * N * M + N - M;
	    
	    double big_tau = 0.2; // big time step of tau leaping
	    const int ratio = int(num_clocks/10); // ratio of big step and small step
	    
	    // From this formula, we can find out that N has to be larger than 10 no matter what
	    // Otherwise, ratio will be zero hence the small_tau would go to infinity
	    double small_tau = big_tau/double(ratio); // small time step
	    
		double* last_update = new double[M * N];
	    
	    random_device rd;
	    
	    mt19937 mt(rd());
	    uniform_real_distribution<double> u(0,1);
	    

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
	    // Keep in mind that time_array is meant to be 1D for convenience.
	    interaction* time_array = new interaction[num_clocks];


	    // Index of clocks (interactions), not the index in the time_array.
		int index = 0;
		
	    for(int n = 0; n <= 2 * N - 2; n++)
	    {	
	    	// If the clock is on the vertical side.
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
			
	        // If the clock is on the horizontal side.
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
	    
	    // Energy profile for testing
	    profile** E_avg_profile = new profile*[N];
	    
	    for(int j = 0; j < N; j++)
	    {
	    	E_avg_profile[j] = new profile[M + 2];
	    }
	    
	    //cout << "Here?" << endl;
	    
	    // Initializing energy profile array
	    for(int i = 0; i < N; i++)
	    {
	    	for(int j = 0; j < M + 2; j++)
	    	{
	    		E_avg_profile[i][j].location = make_pair(i,j);
				// Initialize everything to zero just for convenience
				E_avg_profile[i][j].prev_time = 0.0;
				E_avg_profile[i][j].integral_energy = 0;
			}
		}
		
	    //cout << "index: " << index << endl;
	    
	    int count = 0;
	    
		// each element in the array is the head of a list
	    interaction** clock_time_in_step = new interaction*[ratio + 1];
	    
	    for(int i = 0; i < ratio + 1; i++)
	    {
	        clock_time_in_step[i] = NULL;
	    }
	    
	    gettimeofday(&t1,NULL);
	    
		// big_step_distribute(clock_time_in_step,time_array,N+1,small_tau,ratio,0);

	    
	    int Step = 5;
		   
	    for(int out_n = 0; out_n < Step; out_n++)
	    {
		 	cout << "" << endl;
			cout<< "At Step "<< out_n << endl;
		
	        big_step_distribute(clock_time_in_step, time_array, num_clocks, small_tau, ratio, out_n);
	        
	        for(int in_n = 0; in_n < ratio; in_n++)
	        {
	            if(clock_time_in_step[in_n]!= NULL)
	            {
					cout<< "" << endl;
					cout<< "in_n: " << in_n << endl;
					update(clock_time_in_step, in_n, N, M, small_tau, ratio, out_n, time_array, energy_array, u, mt, &count, E_avg_profile);
					// cout << "Hi1" << endl;
	            }
			// print_v(energy_array,N+2);
	        }
	        clock_time_in_step[ratio] = NULL;
	    }
	    
	    myprofile << endl << "Profile Array: " << "(Size: " << N << "*" << M + 2 << ")" <<endl;
	    myprofile << "===============================================" <<endl;
	    
	    for(int i = 0 ; i < N ; i++) 
	    {
	        myprofile << endl << " " << endl;
			
			for(int j = 1 ; j < M + 1; j++) 
	        {
				myprofile << E_avg_profile[i][j].integral_energy <<" ";
	        }
	        myprofile<<endl;
	    }
		
		cout << " " << endl;
		cout << "Profile Array: " << "(Size: " << N << "*" << M + 2 << ")" <<endl;
	    cout << "===============================================" <<endl;
	    
	    
	    for(int i = 0 ; i < N ; i++) 
	    {
			cout << " " << endl;
			for(int j = 1 ; j < M + 1; j++) 
	        {
				cout << E_avg_profile[i][j].integral_energy <<" ";
	        }
	    }
		
		cout << " " << endl;
		cout << "===============================================" <<endl;
		
	    gettimeofday(&t2, NULL);

		cout << "" << endl;
		cout<< "Done! " <<endl;
	    
	    delete[] energy_array;
	    delete[] E_avg_profile;
	    delete[] time_array;
	    delete[] clock_time_in_step;
	    
	    double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
	                    t2.tv_usec - t1.tv_usec) / 1.e6;
		
		cout << "total CPU time = " << delta <<endl;

	    cout<<" N = "<< N <<endl;
	    cout<<" M = "<< M <<endl;
	    
	    // cout<<"seconds per million event is "<< 1000000*delta/double(count)<<endl;
	    // myfile<<" N = "<<N <<endl;
	    
	    
	    myfile<< 1000000*delta/double(count)<<"  ";
	    myfile.close();

	}

