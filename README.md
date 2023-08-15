# schrodinger_poission_2d_self_consistent_solver
In this project I developed a solver which consecutively solves Poisson and Schrödinger equation to determine the charge as well as potential profile of a double gate (DG) MOSFET. Although double gate MOSFET is much backdated in today’s perspective, it’s a great choice for testing and evaluating the performance of a-self consistent solver.

Double gate MOSFET has a symmetric structure. In this assignment, I have to implement a 1D 
solver for this structure. The only direction of simulation is along the z axis (from upper metal 
gate to lower metal gate in the figure 1.) here. So, the permittivity variation along the z axis 
varies according to this fashion: Metal-Oxide-Channel(body)-Oxide-Metal. To incorporate the 
real-life MOS structure, permittivity of different regions have been chosen from a published 
literature article.

Whatever, our main goal is to solve Poisson equation for the 1st time with a random charge 
density and get the potential profile from it (which is obviously wrong one!). Then it must be
converted to the potential energy I electron volt (eV) which is compatible with Schrödinger
equation. Now, this module (Schrödinger equation) will give wavefunctions as well as E 
profile which is (quantum mechanical equation derived). A new charge density (actually 
carrier concentration) can be derived from the shi and E got from Schrödinger equation
solve. This new charge density (more accurate) give more accurate potential profile hence 
more accurate shi and E. This recursive loop will be running until the solution converges. The 
convergence can be detected by introducing a threshold checkpoint:
 (𝑉𝑛𝑒𝑤 − 𝑉𝑜𝑙𝑑)
2 < 𝑇ℎ𝑟𝑒𝑠ℎ𝑜𝑙𝑑 (1)
The threshold value determines the accuracy and, simulation time. SO, it must be properly 
tunes to trade-off between accuracy and time. The generated results are solved using such a 
threshold value that, it takes 7 iterations to converge the solution.
To determine C-V characteristics, this solver repeats the whole process for different values of 
gate voltage and store the values of corresponding charge against gate voltage. Then the 
derivative of charge to voltage variation can be considered as capacitance and can be 
plotted against gate voltage to depict the C-V characteristics of MOSFET structure for a gate 
voltage regime
