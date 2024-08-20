# Quantum Monte Carlo and Stabilizer States

This is the quantum Monte Carlo code that I wrote to simulate the CNOT model in the paper [*Quantum Monte Carlo and Stabilizer States*](https://arxiv.org/abs/2408.09978). Compile using the command
`g++ -std=c++17  -Ofast -o SSE SSE.cpp`. Starting in line 15, set the number of qubits $N$, the length $L$ of the operator string and the value of the external field $h$ (in units of $J$). 
The number total number of Monte Carlo cycles can be set in line 219 by modifying the `rep_tot` constant.

Please cite the paper if you are using the code.
