# CHP-Simulator
Simulates clifford dynamics with measurements
This code allows for the simulation of the clifford group unitaries which include CNOT gates, phase gates, Hadamard gates and swap gates.  

The code also allows for measurements in the computational basis. 

The entanglement entropy can be calculated with relative efficiency by finding the rank of a binary matrix.

This code is based off results from the paper ["Improved Simulation of Stabilizer Circuits"](https://arxiv.org/abs/quant-ph/0406196) by Scott Aaronson and Daniel Gottesman.

The code file includes a simple example circuit that intializes a state in the |000> state and then creates a GHZ entangled state. The example also determines the entanglement between qubits 1,2 and qubit 3 before and after performing measurements on all three sites.
