//
//
//   The routine(s) in this file are a part of the
//                     C++Lifford
//   suite, developed 2021, and copyrighted
//   to the authors: Cian Reeves and Sagar Vijay
//   at the University of California, Santa Barbara
//
//
//  If you use or modify any part of this routine
//  the header should be kept and unmodified.
//
//  Implements CHP unitaries with measurements and allows for determination of entanglement entropy between regions
//  and is based on the following publication, "Improved Simulation of Stabilizer Circuits", 
//  Scott Aaronson and Daniel Gottesman https://arxiv.org/abs/quant-ph/0406196
//
//

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
