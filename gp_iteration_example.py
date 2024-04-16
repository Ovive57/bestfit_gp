#### Basics ####
import numpy as np

## Program
from gp_setup import sam
from gp_setup import mass_function

from gp_setup import parameter_space

from gp_likelihoods import get_chi2

from gp import gprocesses_model
from gp import fit

## Plots
import matplotlib.pyplot as plt

# Observed values for verifications
obs_acool = 0.8
obs_gammasn = 3.2
obs_properties = sam(obs_acool,obs_gammasn) # returns some useful proerties about the "observed" galaxies, as for example the stellar mass.
obs_function = mass_function(obs_properties[:][0][0]) # returns a function that depends on the observed properties, for example the mass function


# Initial Hammersley
n_init_samples = 10 # number of starting points
dim = 2 # number of dimensions, how many parameters are you varying
limits = [(0.5, 1.5), (0.5, 3.5)] # limits for those parameters, must be len = dim

# Grid to emulate
n_grid = 1000 # grid will be (n_gridxn_grid)
coordinates=parameter_space(n_grid, limits)
coordinates_array = np.array(coordinates)

X_train = np.array(parameter_space(n_init_samples, limits)) # First starting points ("location")
chi2_train = get_chi2(X_train, obs_acool, obs_gammasn) # Actual Chi2 of the starting points

gp_model = gprocesses_model(X_train, chi2_train, dim=dim) # Define the gp_model with the starting points

grid = np.array(parameter_space(n_samples=n_grid, limits=limits))
#plt.plot(grid[:,0], grid[:,1], 'o')
#plt.show()
#exit()

plot = True

for n in range(100):
    chi2, var = fit(grid, gp_model) # Emulate the chi2 and the var for all the points in the grid

    ### POINTS TO ADD:

    ind = np.where(chi2>1e-16) # just save positives chi2

    chi2 = chi2[ind]
    var = var[ind]
    coordinates_array_it = coordinates_array[ind]

    # FIRST NEW POINT TO ADD:
    ind1 = np.argmin(chi2) # Maximum likelihood == minimum chi2 (?)

    # SECOND NEW POINT TO ADD:
    # Weighted sum:
    weight_chi2 = 0.5
    weight_error = 1 - weight_chi2

    weighted_sum = weight_chi2*chi2 + weight_error*var
    ind2 = np.argmax(weighted_sum)


    if plot:
        # PLOT INITIAL POINTS + OBSERVATIONS
        fig,ax = plt.subplots()
        plt.plot(obs_acool, obs_gammasn, '*g', label = 'true point')
        plt.plot(X_train[:,0], X_train[:,1], 'ok', label = 'initial points')

        # PLOT NEW COORDINATES

        plt.plot(coordinates_array_it[ind1,0], coordinates_array_it[ind1,1], 'or', label='new point max likelihood')
        plt.plot(coordinates_array_it[ind2,0], coordinates_array_it[ind2,1], 'ob', label='new point weighted sum')
        plt.legend()
        plt.savefig("plots/new_points_"+str(n)+".jpg")
        plt.close()


    # append the two new points for the next iteration:
    X = np.append(X_train,[coordinates_array_it[ind1]])
    X = np.append(X, [coordinates_array_it[ind2]]).reshape(-1,2)
    X_train = X


    if plot:
        fig,ax = plt.subplots()

        X_plot = coordinates_array_it[:,0]
        Y_plot = coordinates_array_it[:,1]
        Z_plot = chi2
        plt.scatter(X_plot, Y_plot, c=Z_plot, cmap='viridis', vmin=Z_plot.min(), vmax=Z_plot.max())
        ax.scatter(coordinates_array_it[ind1,0], coordinates_array_it[ind1,1], marker='X', c='r', s=100, label='new point max likelihood: ' + str(chi2[ind1]))
        plt.colorbar(label='Chi2 (Z)')
        plt.xlabel('a_cool')
        plt.ylabel('gamma_sn')
        plt.title('Chi2 as Color')
        plt.grid(True)
        plt.legend()
        plt.savefig("plots/chi2_"+str(n)+".jpg")
        plt.close()

        fig,ax = plt.subplots()
        #ind_plot = np.where(var_before!=0)

        Z_plot = weighted_sum

        plt.scatter(X_plot, Y_plot, c=Z_plot, cmap='viridis', vmin=Z_plot.min(), vmax=Z_plot.max())
        ax.scatter(coordinates_array_it[ind2,0], coordinates_array_it[ind2,1], marker='X', c='r', s=100, label='new point weighted sum')

        plt.colorbar(label='var Chi2 (Z)')
        plt.xlabel('a_cool')
        plt.ylabel('gamma_sn')
        plt.title('Weigthed sum Chi2+var as Color')
        plt.grid(True)
        plt.legend()
        plt.savefig("plots/var_chi2_"+str(n)+".jpg")
        #plt.show()
        plt.close()

    # Calculate the chi2 in all the training points: starting points + 2 new points for each iteration
    chi2_train = get_chi2(X_train, obs_acool, obs_gammasn) #! Here is where Shark is run

    # Define a new gp_model with the 2 new points:
    gp_model = gprocesses_model(X_train,chi2_train, dim=dim)


