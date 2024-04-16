#### Basics ####
import numpy as np

from gp_setup import sam
from gp_setup import mass_function

def evaluate_chi2(observed_function, theor_function, error, verbose=False):
    """Calculate the chi2 between the observed and the theoric function

    Args:
        observed_function (array): values for the observed function
        theor_function (array): values for the modeled function
        error(int):(covariance matrix) error for observations and theoretical

    Returns:
        float: chi2
    """
    np.random.seed(42)

    Cn = np.ones((len(observed_function), len(theor_function)))
    error = np.random.rand(len(observed_function), len(theor_function))
    cov = Cn*error
    var = np.diag(cov)
    chi2 = np.sum(((observed_function-theor_function)**2)/(var**2))

    return (1/2)*chi2


def get_chi2(space_array, obs_acool, obs_gammasn, verbose=False):
    """Get chi2 for a set of parameters in a parameter space

    Args:
        space_array (array): parameters in the space of parameters
        obs_acool (float): value for the observed parameter
        obs_gammasn (float): value for the observed parameter
        verbose (bool, optional): _description_. Defaults to False.

    Returns:
        array: array of the chi2 for each value of the parameters
    """
    chi2s = []
    obs_properties = sam(obs_acool, obs_gammasn, verbose=verbose)
    obs_function = mass_function(obs_properties[:][0][0],verbose=verbose)

    for parameters in space_array:
        properties = sam(parameters[0], parameters[1], verbose=verbose)
        function = mass_function(properties[:][0][0],verbose=verbose) # Taking only the first one as it is the mass for example
        chi2 = evaluate_chi2(obs_function, function,0.2,verbose=verbose)
        chi2s.append(chi2)

    return np.array(chi2s).reshape(-1,1)