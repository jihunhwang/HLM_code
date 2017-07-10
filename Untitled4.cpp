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

struct interaction
{
    
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
    
};


/*
Create the time array (an array of clocks (interactions))
*/
interaction* make_time_array(const int N, const int M)
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
	interaction return_array[2 * N * M + N - M];
	
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
	
	/*
	for(int i = 0; i < 2N - 1; i++)
	{
		if(i % 2 == 0)
		{
			return_array[]
		}
	}
	
	for(int i = 0; i < num_clocks; i++)
	{
		return_array[i].index = i;
	}
	*/
}

int main()
{
	interaction* test_array = make_time_array(5, 5);
	
	
	for(int a = 0; a < 50; a++)
	{
		cout << test_array[a].index << endl;
	}
	
	
	/*
	for(int a = 0; a < 50; a++)
	{
		cout << test_array[a].location.first << endl;
	}*/
	
	
	/*
	for(int a = 0; a < 50; a++)
	{
		cout << test_array[a].ishorizontal << endl;
	}*/
	
}


