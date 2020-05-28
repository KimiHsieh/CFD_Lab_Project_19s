# CFD Lab Final Project #
This repository contains *Master-Praktikum: Scientific Computing - Computational Fluid Dynamics* final project for Technische Universität München (TUM) [IN2186, IN2106] course in summer 2019.

The coursework description is available [here](https://campus.tum.de/tumonline/wbLv.wbShowLVDetail?pStpSpNr=950461641&pSpracheNr=2&pMUISuche=FALSE). The full project report is available [here](https://drive.google.com/file/d/1svNHnDHDr8D6bwV7smY62-gkLlsaducp/view?usp=sharing)

# Introduction #
In this project we implemented diﬀerent time-stepping methods for problems that apply Navier-Stokes Equations, considering the geometry of the domain. The scope is to analyse if using a more precise integration algorithm is necessary and achievable with a not very much increased cost. The numerical integration methods implemented in the current version of the project are:
## Fixed time-step ##
* Explicit Euler
* Heun method (Runge-Kutta 2nd order)
* Runge-Kutta 4th order (classical Runge-Kutta)

## Adaptive time-step ##

* Bogacki–Shampine (Runge-Kutta 3/4) (MATLAB: ode23)

* Runge–Kutta–Fehlberg (Runge-Kutta 4/5) (MATLAB: ode45)

This repository contains

- a makefile
- the headers
- the files with the respective method stubs
- a folder with parameter files (dat)
- a folder with geometry files (geometry)


1. Running a simulation

Use the format:	 ./sim [problem_name],
where [problem_name] is one of our defined problems:
a) cavity100 			 stationary problem: Driven cavity (analyzed in WS1)
b) karman_vortex 		 non-stationary problem: The Karman Vortex Street (analyzed in WS2)
c) channel-bfs 			 stationary problem: Flow over a Step (analyzed in WS2)
or your own data file put in the 'dat' folder.
Moreover, take care that the geometry files should be put in the 'geometry' folder, as pgm files.
e.g.: ./sim karman_vortex
Important! Do not modify the structure of the .dat files, because of the hard-coding used for compatibility in reading the parameters.

1. General remarks

This version of the project implements the Navier-Stokes equations in 2D for arbitrary geometries. The coupled variables are then the two components of the velocity and the pressure (U,V,P).
The program takes two input files, data (.dat) and geometry (.pgm), that contain the problem parameters and the computational domain. They are stored respectively in the folders 'dat' and 'geometry'.
The outputs are vtk files used for visualization in Paraview and put in the ./paraviewFiles/[problem name] folder.
To run the code, one needs to use the executable sim and specify the particular problem (./sim [problem name]). If an attempt is made to run 'sim' without specification of the particular problem, a message shows up telling the user what the names of the available problems are. No recompilation is needed to simulate a different problem.

2. Implemented time-stepping methods
We implemented 5 time-stepping methods which are selected from the configuration file by using the parameter method
Fixed time-step methods:
- Explicit Euler (former method)                          (method = 0)
- Heun method (Runge-Kutta 2nd order)                     (method = 1)
- Runge-Kutta 4th order (classical Runge-Kutta)           (method = 2)
Adaptive time-step methods:
- Bogacki–Shampine (Runge-Kutta 3) (MATLAB: ode23)        (method = 3)
- Runge–Kutta–Fehlberg (Runge-Kutta 4/5)  (MATLAB: ode45) (method = 4)
For the adaptive time-step methods we have an additional parameter - toler (tolerance) - used to decide when to adapt the time-step. We set it to 10%, i.e. if the error of estimator is greater than 10% we adapt the time-step.


Contributed by Yi-Han Hsieh, Teodor Rotaru, Mario Ponce Mart´ınez