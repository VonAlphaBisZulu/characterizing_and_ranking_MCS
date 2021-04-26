# Characterizing and Ranking Computed Metabolic Engineering Strategies
The files in this directory contain a script and the necessary functions to reproduce the examples from:

"Characterizing and Ranking Computed Metabolic Engineering Strategies" (submitted, Dec 2018)
Philipp Schneider, Steffen Klamt
Max Planck Institute for Dynamics of Complex Technical Systems, Sandtorstraße 1, 39106 Magdeburg, Germany

-------------

File list:

rankExample.m
characterizeIS.m
iJO1366-1-4-Butanediol.mat
iJO1366-L-Methionine.mat
cell2csv.m

-------------

Requirements:

CellNetAnalyer	(2018.1 or above)
(you may download it from: http://www2.mpi-magdeburg.mpg.de/projects/cna/cna.html)

MATLAB 		(2013 or above)

IBM CPLEX 		(12.6.3 or above)
(see CellNetAnalyzer manual for instructions regarding the interplay von CPLEX
and CellNetAnalyzer)

-------------

How to:

 0.  Unzip the file "RankingScripts.zip". (A directory "RankingScripts" is generated).
 1.  Start MATLAB from the CellNetAnalyzer directory.
 2.  Start CellNetAnalyer (within MATLAB) with "startcna"
 3.  Add the path to the "RankingScripts" directory via addpath('....../RankingScripts').
(4.) Specify the example to be computed via the 'filename' variable in 'rankExample.m'.
 5.  Run 'rankExample.m' (from the CellNetAnalyzer main directoy as current working folder in MATLAB).

The "rankExample" script uses only the first four MCS to demonstrate the
ranking procedure. If you want to rank all MCS (and fully reproduce the 
results fom the paper) please outcomment the respective (and marked) line! 

-------------

Content:

1. **rankExample.m**

Running this MATLAB script (after starting CellNetAnalyer) will reproduce the results presented in the work mentioned above. It loads the required models, minimal cut sets (MCS) and variables from a '.mat' file and prepares all data to fit the format required in the actual characterization and ranking function 'characterizeIS'.

2. **characterizeIS.m**

This MATLAB function characterizes a given set of intervention strategies (IS) among up to 10 different criteria. If specified, this function also classifies the given set of ISs into equivalence classes and provides a scoring and ranking of the candidates. A detailed description is provided in the function's documentation.

3. **iJO1366-1-4-Butanediol.mat**

Is a container for several variables:
cnap:		the E. coli ijO1366 model that contains additional reactions to produce 1,4-Butanediol.
mcs:		a set of minimal reaction cut sets that enforce strong growth coupling
D,d:		vectors that specify which flux distributions must remain feasible 	(feasible:   D*r <= d)
T,t:		vectors that specify which flux distributions must remain infeasible 	(infeasible: T*r <= t)
coreReacs:	A vector containing the reaction indices used to generate the Core model.

4. **iJO1366-L-Methionine.mat**

Is a container for several variables:
cnap:		the E. coli ijO1366 model that contains additional reactions to produce L-Methionine.
mcs:		a set of minimal reaction cut sets that enforce strong growth coupling
D,d:		vectors that specify which flux distributions must remain feasible 	(feasible:   D*r <= d)
T,t:		vectors that specify which flux distributions must remain infeasible 	(infeasible: T*r <= t)
coreReacs:	A vector containing the reaction indices used to generate the Core model.

5. **cell2csv.m**

A function that saves a cell matrix to the csv format. the function was retrieved from the MATLAB file exchange repository. Authors and function history are enclosed in the file.
