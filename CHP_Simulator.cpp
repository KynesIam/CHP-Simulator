#include <iostream>
#include <algorithm> 
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>

/*
    Implements CHP unitaries with measurements and allows for determination of entanglement entropy between regions.


    Reference:
        "Improved Simulation of Stabilizer Circuits"
        Scott Aaronson and Daniel Gottesman
        https://arxiv.org/abs/quant-ph/0406196
*/
std::random_device rd;
std::mt19937 gen(rd());

 
    
    
void initialise_state(std::vector<std::vector<int>> &state)
{
        /*Initializes state to be in the computational basis
         
         
          Arguments:
                    State: The tableau representation of the state
        */
        {
            for(int i = 0; i<state[0].size(); ++i)
            {
                    state[i][i] =  1;
            }
        }
}

        
class Chp_simulator
{
    

    public:
        int L;
        std::vector<std::vector<int>> state;
        std::vector<int> r;
        Chp_simulator(int number_of_qubits, std::vector<std::vector<int>> initial_state, std::vector<int> initial_r)//Initialisation of state and stabiliser signs
        {
            L = number_of_qubits;
            state = initial_state;
            r = initial_r;
        }
        
        
        
        void hadamard(int site)
        {
            /*Performs hadamard gate on the site "site"
            
            
            Arguments: 
                    site: site to which hadamard gate is applied
            */
            
            for (int a = 0;a<2*L; a++)
            {
                r[a] ^= state[a][site]*state[a][site+L];
                std::swap(state[a][site],state[a][site+L]);
            }
        }

        void cnot(int control, int target)
        {
            /*Performs cnot gate between 2 sites
             
             Arguments:
                        control: the control site for the CNOT gate
                        target:  the target site for the CNOT gate
             */
            for (int a = 0;a<2*L; a++)
            {
                r[a] ^= (state[a][control]*state[a][target+L])*(state[a][target]^state[a][control+L]^1);
                state[a][target]^= state[a][control];
                state[a][control+L]^= state[a][target+L];
            }
        }

        void phase(int site)
        {
            /*Perfroms a phase gate(S gate) to a site
             
             
             Arguments:
                        site: the site to which the gate is applied
             */
            for (int a = 0;a<2*L; a++)
            {
                r[a] ^= state[a][site]*state[a][site+L];
                state[a][site+L] ^= state[a][site];
            }
        }

        void swap(int site1, int site2)
        {
            /*Perfroms a swap between two sites
             
             
             Arguments: 
                        site1: one of the sites to be swapped
                        site2: the second site to be swapped
             */
            for (int a = 0;a<2*L; a++)
            {
                std::swap(state[a][site1],state[a][site2]);
                std::swap(state[a][site1+L],state[a][site2+L]);
            }
        }

        int measure_z(int site)
        {
            /*Performs a measurement in the Z-basis
             
             
             Arguments:
                        site:  The site to be measured
                        
            Returns:
                    Result of measurement on the chosen site, either 0 implying |0> state or 1 implying |1> state.
             */
            
            std::uniform_int_distribution<> outcome(0, 1);
            
            for(int p = L; p<2*L; p++)
            {
                if(state[p][site] == 1)//State has X-stabiliser and therfore measurement outcome will be probabilistic
                {
                    //Following updates the remaining stabilizers ensuring state[p][site] is the only stabilizer with an X at the chosen site
                    for(int i = 0; i<2*L;i++)
                    {
                        if(state[i][site] == 1 and i != p)
                        {
                            row_sum(i,p);
                        }
                    }
                    
                    for(int i = 0; i<2*L;i++)//Swap the p and (p-L)th row and set the p-th row to 0
                    {
                        state[p-L][i] = state[p][i];
                        state[p][i] = 0;
                    }
                    
                    state[p][site+L] = 1;//initialise this site to be in Z-state
                    r[p] = int((1+pow(-1,outcome(gen)))/2);//initialise it to be |0> or |1> state with equal probability
                    return (r[p]); 
                    
                }
            }
            
            //Below deals with case where the measurement outcome is already determined
            for(int i = 0; i<2*L;++i)
            {
                state[2*L][i] = 0;
            }
            r[2*L] = 0;
            for(int i = 0; i<L;i++)
            {
                if(state[i][site] == 1)
                {
                    row_sum(2*L,i+L);
                }
            }
            return r[2*L];
                
        
        }
        
        int g(int x1, int z1, int x2, int z2)
        {  
            /*Used to calculate n in (i)^n when multiplying two pauli matices.
             
              e.g X*Y = iZ so n = 1
              
              Arguments:
                        x1: Determines whether first Pauli matrix contains an X
                        z1: Determines whether first Pauli matrix contains a Z
                        x2: Determines whether second Pauli matrix contains an X
                        z2: Determines whether second Pauli matrix contains a Z
            */
            
            if(x1 == 1 and z1 == 1) //Y gate.
            {
                //No phase for YI = Y
                //-1 phase for YX = -iZ
                //No phase for YY = I
                //+1 phase for YZ = +iX
                return z2-x2;
            }
            
            else if(x1 == 1 and z1 == 0)//X gate.
            {
                //No phase for XI = X
                //No phase for XX = I
                //+1 phase for XY = iZ
                //-1 phase for XZ = -iY
                return z2*(2*x2-1);
            }
            
            else if(x1 == 0 and z1 == 1)//Z gate.
            {
                //No phase for ZI = Z
                //+1 phase for ZX = -iY
                //-1 phase for ZY = iX
                //No phase for ZZ = I
            
                return x2*(1-2*z2);
            }
            
            //x1 = z1 = 0
            return 0;
                
        
        }

        void row_sum(int h, int j)
        {
            /*Multiplies two pauli strings together
             
             
             Arguments:
                        h: row corresponding to first pauli string in "state" tableau
                        j: row corresponding to second pauli string in "state" tableau
             */
            int m = 2*(r[h] + r[j]);
            for(int k = 0; k<L; k++)
            {
                m+= g(state[j][k],state[j][k+L],state[h][k],state[h][k+L]);
                state[h][k] ^= state[j][k];
                state[h][k+L] ^= state[j][k+L];
            }
            r[h] = int((m%4)/2);
        }

        void print_tableau()
        {
            for(int i = 0; i<2*L+1; ++i)
            {
                for(int j = 0; j<2*L; ++j)
                {
                    std::cout<<state[i][j]<<" ";
                }
                std::cout<<"\n";
            }
        }
        
        void print_stabilizers()
        {
            //Prints out the explicit stabilizers for the state including the sign
            for(int i = L; i<2*L; ++i)
            {
                std::string stabilizer;
                if(r[i] == 0)
                {
                    stabilizer+="+";
                }
                else if(r[i] == 1)
                {
                    stabilizer+="-";
                }
                
                for(int j = 0;j<L;++j)
                {
                    if(state[i][j] == 1 and state[i][j+L] == 1)
                    {
                        stabilizer+="Y";
                    }
                    else if(state[i][j] == 0 and state[i][j+L] == 1)
                    {
                        stabilizer+="Z";
                    }
                    else if(state[i][j] == 1 and state[i][j+L] == 0)
                    {
                        stabilizer+="X";
                    }
                    else if(state[i][j] == 0 and state[i][j+L] == 0)
                    {
                        stabilizer+="_";
                    }
                }
                std::cout<<i-L<<":"<<stabilizer<<"\n";
            }
            std::cout<<"\n";
        }
        
        std::vector<std::vector<int>> get_A_stabilizers(int region_start,int region_end)
        {
            /*Gets the stabilizers for a contiguous region of size |region_start - region_end|
                
                
              Arguments:
                        region_start: starting point of region for which stabilizers are found
                        region_end:  finishing point of region for which stabilizers are found
            */
            
            std::vector<std::vector<int>> region_A(2*(region_end-region_start), std::vector<int>(L));
            for(int i = region_start; i<region_end; ++i)
            {
                for(int j = L; j<2*L; ++j)
                {
                    region_A[2*i][j-L] = state[j][i];
                    region_A[2*i+1][j-L] = state[j][i+L];
                }
            }
            return region_A;
        }
};

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


