#### Basics ####
import numpy as np
import time

#GPY
import GPy

# Likelihood minimisation
from scipy.stats import multivariate_normal

def gprocesses_model(X, y, dim, kernel_type = "matern52", num_restarts = 5, verbose=False):
    """Define the gaussian processes model given a kernel and a training set.

    Args:
        X (array): array of "coordinates". Te values of the parameters
        y (array): array of the likelihoods (chi2) associate a these parameters
        dim (int): dimension of our problem. Number of parameters that we are optimising.
        kernel_type (str, optional): name of the kernel to use. Defaults to "matern52".
        num_restarts(int, optional): number of restarts to avoid a minimum. Defaults to 5.
        verbose (bool, optional): if information is needed. Defaults to False.

    Returns:
        class 'GPy.models.gp_regression.GPRegression'
        """

    if kernel_type == "matern52":
        # Define Matérn kernel with length scale between 10e-3 and 10
        kernel = GPy.kern.Matern52(input_dim=dim, variance=1.0, lengthscale=1.0)
        # Set constraints on the length scale parameter:
        # Limits following Rocher+23:
        kernel.lengthscale.constrain_bounded(1e-3, 10)
    else:
        print(f"WARNING: Kernel {kernel_type} not supported yet.")
        exit()

    if verbose:
        # Print kernel summary
        print("\nKERNEL\n\n", kernel)


    gp_model = GPy.models.GPRegression(X,y, kernel=kernel)
    if verbose:
        print("\n\nGP MODEL\n ", gp_model)

    #! This does not work:
    # Specify the likelihood function (Gaussian likelihood)
    #GP_likelihood = GPy.likelihoods.Gaussian()

    # Set the likelihood function for the GP model
    #gp_model.likelihood = GP_likelihood

    # Specifying the optimization method (L-BFGS-B) and number of restarts
    # L-BFGS-B stands for Limited-memory Broyden–Fletcher–Goldfarb–Shanno with Bound constraints.
    gp_model.optimize('lbfgsb', messages=verbose)
    gp_model.optimize_restarts(num_restarts=num_restarts, messages=verbose)
    if verbose:
        # Print gp_model summary
        print("\n\nGP OPTIMISED MODEL\n ",gp_model)
    return gp_model


def fit(grid, gp_model, verbose=False):
    """Fit a multivariate Gaussian to the emulated likelihood and estimate the mean and covariance.

    Args:
        grid (array): points at which to make the prediction
        gp_model ('GPy.models.gp_regression.GPRegression'): gaussian processes model
        verbose (bool, optional): _description_. Defaults to False.

    Returns:
        mean, var: the predicted chi2 and variance in each of the points in the grid

    """

    # full_cov = False means, we get only the diagonal of the cov matrix, i.e. the variance
    # include_likelihood (bool) – Whether or not to add likelihood noise to the predicted underlying latent function f.
    res = gp_model.predict(grid, full_cov=False, include_likelihood=False)
    mean, var = (res[0].T)[0], (res[1].T)[0]

    for i, m in enumerate(mean):
        if np.abs(m)<1e-16:
            mean[i]=0
    ind = np.where(mean!=0)

    #https://www.statlect.com/fundamentals-of-statistics/normal-distribution-maximum-likelihood
    # In Rocher et al (2023): "As often assumed in Gaussian processes, the mean value of the joint Gaussian distribution is assumed to be 0.""
    # so, if the likelihood is proportional to exp(-1/2 ((x-nu)²/var)), we assume x = mean (emulated) and nu = 0
    # I am not sure about this explanation, maybe here is the problem

    # following the same structure as in gp_likelihoods for the chi2: chi2 = np.sum(((observed_function-theor_function)**2)/(var**2))
    chi2 = np.exp((1/2)*((mean**2)/var**2))


    #return mean, var #! returns mean and not chi2
    return chi2, var