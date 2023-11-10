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
#include <cmath>



int entanglement_entropy(std::vector< std::vector<int>> stabilizers)
{
    /*finds the rank of the input array.  The entanglement entropy of a region A is then given by rank(A) - |A|, where A is an array of the stabilizers on region A
    
    Arguments: 
            stabilizers: The stabilizers of the desired state
    
    */        
    int n = stabilizers[0].size();
    int rank = 0;
    
    for (int col = 0; col< n;col++)
    {
            int j = 0;
            std::vector<int> rows(0);
            while(j<stabilizers.size())
            {
                if (stabilizers[j][col] == 1)
                {
                    rows.insert(rows.end(), j);
                }
                j+=1;
            }
            if(rows.size()>=1)
            {
                for(int c = 1; c<rows.size(); c++)
                {
                    for(int k = 0; k<n;k++)
                    {
                        stabilizers[rows[c]][k] ^= stabilizers[rows[0]][k];
                        
                    }
                    
                }
                stabilizers.erase(stabilizers.begin()+rows[0]);
                rank+=1;
            }
    }
    
    for(int i = 0; i<stabilizers.size(); i++)
    {
        int sum = 0;
        for(int j = 0; j<stabilizers[i].size();j++)
        {
            sum += stabilizers[i][j];
            
        }
        if(sum>0)
        {
            rank+=1;
        }
    }
    return rank;
}
