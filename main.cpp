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
#include <algorithm> 
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include "./headers/circuit.h"
#include "./headers/entanglement_entropy.h"
#include "./headers/state_initialization.h"

int main()
{
    //Simple example circuit: Creates GHZ state then measures all 3 sites and calculates the entanglement entropy between site 12 and site 3.
    //0:-----------X---M---
    //             |
    //1:-------X---|---M---
    //         |   |
    //2:---H---#---#---M---
    
    int L = 3;
    std::vector<std::vector<int>> initial_state(2*L+1, std::vector<int>(2*L));
    std::vector<int> initial_r(2*L+1);
    initialise_state(initial_state);
    
    Chp_simulator circuit(L, initial_state,initial_r);
    
    circuit.hadamard(2);
    circuit.cnot(2,1);
    circuit.cnot(2,0);
    
    circuit.print_stabilizers();
    std::cout<<"Entanglement entropy before measurement: "<<entanglement_entropy(circuit.get_A_stabilizers(0,2)) - 2<<"\n";
    std::cout<<std::endl;
    
    circuit.measure_z(0);
    circuit.measure_z(1);
    circuit.measure_z(2);
    
    circuit.print_stabilizers();
    std::cout<<"Entanglement entropy after measurement: "<<entanglement_entropy(circuit.get_A_stabilizers(0,2)) - 2<<"\n";
    std::cout<<std::endl;

    return 0;
}


