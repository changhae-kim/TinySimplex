"SimplexUtil.py" defines functions and variables which are used in the other scripts. Most of the codes are concerned with I/O and file management.

"SimplexInit.py" initializes the NMS method by generating the vertices. As the input, it takes a Gaussian COM file, which has been edited to mark the variable Z-Matrix coordinates with exclamations (!). See "DCzTRZ.com" for example. As the output, it creates Q-Chem IN files for the initial vertices and a "simplex.log" file which stores information about the vertices.

"SimplexAuto.py" contains the main workflow of the NMS method. At the start of each iteration, it reads the "simplex.log" file and the Q-Chem OUT files for the vertices. Then, it performs reflection, expansion, and contraction (REC) to generate new candidate vertices and runs Q-Chem. At the end of the iteration, it writes the new vertices to the "simplex.log" file.

Some parts of the scripts assume that they can interact with the OS in specific ways, so they might need to be modified to work on other computers.
