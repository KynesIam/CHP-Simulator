#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>


void initialise_state(std::vector<std::vector<int>> &state)
{
    /*
    Initializes state to be in the computational basis
         
    Arguments:
        State: The tableau representation of the state
    */
        for(int i = 0; i<state[0].size(); ++i){
                state[i][i] =  1;
        }
}
