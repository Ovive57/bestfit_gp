# bestfit_gp: How to use

## gp_setup.py

Module that prepares the parameter space. It also has a function that works as a semi-analytical model and another function that works as a mass function. These last two functions will need to be replaced by the actual SAM and the actual function or functions that will be used.

## gp_likelihoods.py

Module that calculates the chi-square of the values of the function used given by the parameters. Normally it is only used for evaluate the function of the initial parameters. Once convergence is implemented (if applicable), this module may be used to recalculate additional chi-squared values for the newly added points.

## gp.py

Module that builds the GP model and fits a multivariate Gaussian to the emulated likelihood. This function returns the mean and the variance. There is a link in line 89 that explains what to do with it to calculate the chi-squared. However, **the code is not working well**. The values it calculates don't make sense. It returns an array where the first value is significantly different from the rest of the values. This causes the algorithm to remain stuck in this value indefinitely. Suggestion: try excluding the first value of the array to see if the algorithm works without it.


## gp_iteration_example.py

Script that performs 100 iterations of the algorithm and generates plots if necessary. There are variables that can be adjusted to conduct additional tests.


# Aditional things:

- In plot_functions/plot_functions.py, you can find some functions that calculate and plot the functions that could be used in the algorithm. Please note that the paths need to be changed and that you need additional data: the SAM data and the observed data.
- To run galform: > cd galform/build
                  > ./galform2 /path/to/output/folder/ /path/to/config/file
- To run Shark: > cd shark/build
                > ./shark /path/to/config/file
