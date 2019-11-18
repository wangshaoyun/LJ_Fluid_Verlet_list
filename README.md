![Schmetical Diagram](https://github.com/wangshaoyun/LJ_Fluid_Verlet_list/blob/master/gr.jpg "Simulation System")
# Monte Carlo Simulation of Lennard Jones Fluid by Verlet Lists Method
## About the program
1. This a program that simulate the equation of state and radial distribution function of LJ fluid by the combination of Verlet lists method and cell list method [1]. 
2. This program is used to very the result of pure cell lists method forwarded by me.
>[1] Frenkel D, Klein M, Parrrinello M. "Understanding Molecular Simulation: From Algorithm to Applications". 2nd ed. Singapore: Elsevier Pte Ltd, 2002.  
>[2] S. Y. Wang, C. H. Tong. Cell-lists Method for Monte Carlo Simulation, to be submitted. 
 
## About Parameter Files 
+ energy_data.txt: It is a parameter file which includes Bjerrum length, FENE parameters, bending stiffiness and torsion stiffness.
+ system_data.txt: It is a parameter file which includes system, time, brushes parameter.  

## Compiling files
1. Determine the parameters such as box size, brush chians, time step and so on.
2. Set the parameters in energy_data.txt and system_data.txt
3. Open terminal and enter the current file content
4. To compile the files:
```
$ make
```
  
5. To execute the files:
```
$ ./main &
```





