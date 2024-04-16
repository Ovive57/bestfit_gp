#### Basics ####
import numpy as np

#### Parameter space ####
from skopt.space import Space
from skopt.sampler import Hammersly



def parameter_space(n_samples, limits, skip=-1, space_type="Hammersley", verbose=False):
    """Generate a parameter space with Hammersly distribution

    Args:
        n_samples (int): number of points in the parameter space
        limits (array): bounds for each parameter in the parameter space
        skip(int): changes the order of one of the parameters, but not generates a new one.
    Returns:
        tupple: parameter space
    """
    space = Space(limits)
    if space_type=="Hammersley":
        hammersly = Hammersly(min_skip=skip, max_skip=skip)
        ps = hammersly.generate(space.dimensions, n_samples, random_state=None)
    return ps


def sam(a,b,nvalues=4, ndata=50, verbose=False):
    """model that will be the semianalitical model

    Args:
        a (array): first parameter
        b (array): second parameter
        nvalues (int): number of properties that I want. Defaults to 4. Number of columns that I am going to have. Only purpose: realism: for example, mass, luminosity, SFR, R50.
        ndata (int): number of data for each property, "size of the column". 50 values, 50 galaxies.

    Returns:
        array: the modelisation
    """
    np.random.seed(42)

    # Check if inputs are floats or arrays
    if isinstance(a, (float, int)) and isinstance(b, (float, int)):
        # If inputs are floats, convert them to arrays
        a = np.full(1, a)
        b = np.full(1, b)
    mod = []

    for ai, bi in zip(a, b):
        inner_arrays = []
        for _ in range(nvalues):
            random_array = np.random.uniform(ai, bi, size=(ndata,))
            inner_arrays.append(random_array)
        mod.append(inner_arrays)

    return mod

def mass_function(prop, verbose=False):
    """function where we will calculate the likelihood

    Args:
        prop (array): first property
    """

    function = np.sin(prop) #+ np.random.normal(0,0.1, size=prop.shape)

    return function