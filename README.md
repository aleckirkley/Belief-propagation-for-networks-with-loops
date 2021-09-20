# Belief propagation for networks with loops
C++ code for computing average magnetization (M), specific heat (C), and entropy (S) for the zero-field Ising model using the neighborhood approximation in 'Belief propagation for networks with loops'. Ipython notebook with Python wrapping provided for easy use. Example in notebook given for power grid network used in Fig. 3 of paper. 

Prior to running Python wrapper code, cd to loopy_BP directory and compile C++ code with 'make'.

Functions for computation of most physical quantities defined in 'include/neighborhoods.h'. These can be modified to adapt algorithm for use with other graphical models beyond the Ising model.

If you use this algorithm, please cite:

A. Kirkley, G. T. Cantwell, and M. E. J. Newman, Belief propagation for networks with loops. Science Advances 7, eabf1211 (2021)

 
