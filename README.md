# Quantum Monte Carlo and Stabilizer States

This is the quantum Monte Carlo code that I wrote to simulate the models in the paper [*Quantum Monte Carlo and Stabilizer States*](https://arxiv.org/abs/2408.09978). SSE.cpp is the Quantum Monte Carlo code to simulate the CNOT model. SSE_TFI.cpp is the Monte-Carlo code for Transverse-Field Ising (TFI) model.
Compile using the command `g++ -std=c++17  -Ofast -o SSE SSE.cpp`. Starting in line 15, set the number of qubits $N$, the length $L$ of the operator string and the value of the external field $h$ (in units of $J$). 
The total number of Monte Carlo cycles can be set in line 219 by modifying the `rep_tot` constant.

The file TFI_ED_Notebook.ipynb is a Jupyter Notebook containing the exact diagonalization results of the TFI model. CNOT_ED.py is the exact diagonalization code for the CNOT model.

Please cite the paper if you are using the code.
