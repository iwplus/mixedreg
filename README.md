# mixedreg
Mixed non-parametric regression between Kernel and Spline Truncated with Particle Swarm Optimization (PSO) for optimum bandwidths and knots search

Use **indeks.csv** file to specify which independent variables follow _Kernel Regression_ or _Spline Truncated_. The 'x' column for kernel and 't' column for spline. Put 1 if it follows the specific regression (either kernel or spline), otherwise put 0. 

_Note_: The number of non-empty rows in both columns must be equal in order R can read the indeks.csv file correctly
