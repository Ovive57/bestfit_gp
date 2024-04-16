##### "Halo mass function"(lhm) and "Stellar mass function"(lm) from a .hdf5 file. #####

### FILES ###
file_shark = '../shark_output/UNIT-PNG/subf+subl+dhalos/128/multiple_batches/galaxies.hdf5' # redshift 0
file_galform = '../galform_output/galaxies.hdf5'

### COSMOLOGY ###
#cosmo_shark = h5ls galaxies.hdf5/cosmology

# Histograms

import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import style

plt.style.use(style.style1)

#f = h5py.File(file_shark,'r')
#run_info = f['run_info']
#volume = run_info['lbox'][()]**3


def plot_lm(model, dm,z=[0,1,2]):
    """Plot the lm for the different models and compare them all with the observations in the same plot for different redshifts.

    Args:
        model (array of str): name of the models
        dm (float): bin of mass for the model's histograms
        z (array of float or int): the different redshifts. Defaults to [0,1,2].
    """
    h = 0.6774
    for iz in z:
        fig,ax=plt.subplots()
        for mod in model:
            if mod == 'shark':
                if iz == 0:
                    filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/128/multiple_batches/galaxies.hdf5'
                if iz == 1:
                    filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/97/multiple_batches/galaxies.hdf5'
                if iz == 2:
                    filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/78/multiple_batches/galaxies.hdf5'

            if mod == 'galform':
                filename = '../galform_output/galaxies.hdf5'

            with h5py.File(filename,"r") as f:
                if mod == 'shark':
                    volume = 1.56250E+07 #f['run_info/lbox'][()]**3 # [Mpc/h]³ #! es 1e9 el 1.56250E+07 es el volumen de 1 subvolumen, 1e9/64=1.56250E+07
                    #mres = f['run_info/particle_mass'][()] # [Msun/h] #! Es 0, mirar por qué
                    mres = 9.97e9 # [Msun/h]
                    mvir = f['galaxies/mvir_hosthalo'][:] # Dark matter mass of the host halo in which this galaxy resides [Msun/h]
                    mstars = f['galaxies/mstars_bulge'][:] + f['galaxies/mstars_disk'][:] #[Msun/h]
                    #print(len(mstars))

                if mod == 'galform':
                    volume = 1.56250E+07 #1e09 # [Mpc/h]³ # more used-parameters
                    mres = 9.97e9 # [Msun/h]

                    mstars = f['Output001/mstars_bulge'][:] + f['Output001/mstars_disk'][:] #[Msun/h]
                    #print(len(mstars)) #! Hay muy pocas, porque tengo galform corrido en 1 solo subvolumen

                lm = np.log10(mstars)
                #print(np.shape(lm))
                ind = np.where(mvir>(mres*20)) # All of them are bigger
                #print(mres*20)
                lm = lm[ind]
                #print(mvir.min(), mvir.max())
                #print(np.shape(lm))
                #exit()
                #mmin = np.log10(mres*10)
                mmin = lm.min() #mres = galaxies.hdf5/run_info/particle_mass
                mmax = lm.max()
                edges = np.array(np.arange(mmin, mmax+dm, dm)) #from mmin to mmax with a dex bin
                hist = edges[1:]-0.5*dm
                nm = np.zeros(shape=len(hist))

                iline=0

                for mass in lm: # in log

                    for ie, edge in enumerate(edges[:-1]):
                        if mass >= edge and mass<edges[ie+1]:
                            nm[ie] = nm[ie] + 1

                    iline+=1

                #print(sum(nm))

                f.close

            yhist = np.log10(nm) - np.log10(dm) - np.log10(volume)

            ax.plot(hist, yhist, label = mod)
            plt.title('Stellar Mass Function z = ' + str(iz))

            plt.ylabel('log$_{10}$ (dn/dlog (M$_{*}$)/h$^{3}$ Mpc$^{-3}$)')
            plt.xlabel('log$_{10}$(M$_{*}$/h$^{-1}$ M$_{\\odot}$)')


        ### OBSERVATIONS ###
        # henriques: [h^-2Msun] If one wants to plot quantities without h in the units, then both observations and simulations should be divided by the h^2 of the simulation.
        # baldry: [Msun]
        if iz==0:
            # Observations henriques
            fileobs='../Obs_Data/smf/henriques_2014_z0_cha.txt'
            mass_low = np.loadtxt(fileobs, skiprows=5, usecols=(0), unpack=True) - np.log10(h)
            mass_high = np.loadtxt(fileobs, skiprows=5, usecols=(1), unpack=True) - np.log10(h)

            dm_obs = mass_high[1]-mass_low[1] # 0.25
            #print(dm_obs) # = 0.25
            hist_obs = mass_high - 0.5 * dm_obs
            yhist_obs = np.log10(np.loadtxt(fileobs, skiprows=5, usecols=(2), unpack=True))#*dm)
            yhist_err = np.abs(np.log10(np.loadtxt(fileobs, skiprows=5, usecols=(3), unpack=True)))#*dm)

            plt.plot(hist_obs, yhist_obs, 'o', label = 'henriques_2014_z' + str(iz))
            #plt.errorbar(hist_obs, yhist_obs, yhist_err,fmt='none', elinewidth=None, ecolor = 'orange')

            # Observations Baldry
            fileobs = '../Obs_Data/smf/baldry_2012_z0_cha.txt'
            mass = np.loadtxt(fileobs,skiprows = 3, usecols = (0), unpack = True) + np.log10(h)
            phi = np.loadtxt(fileobs,skiprows=3, usecols = (1), unpack = True)
            yhist_obs = np.log10(phi)-3*np.log10(h)-3
            yhist_err = np.abs(np.log10(np.loadtxt(fileobs,skiprows=5, usecols = (3), unpack = True))-3*np.log10(h)-3)

            plt.plot(mass, yhist_obs, 'o', label = 'baldry_2012_z' + str(iz))
            #plt.errorbar(mass, yhist_obs, yhist_err,fmt='none', elinewidth=None, ecolor = 'limegreen')


            # Observations Moustakas
            fileobs = '../Obs_Data/smf/moustakas_z0.01_z0.20.smf'
            mass_low = np.loadtxt(fileobs, skiprows=5, usecols=(0), unpack=True) - np.log10(h)
            mass_high = np.loadtxt(fileobs, skiprows=5, usecols=(1), unpack=True) - np.log10(h)

            dm_obs = mass_high[1]-mass_low[1] # 0.25
            #print(dm_obs) # = 0.25
            hist_obs = mass_high - 0.5 * dm_obs
            yhist_obs = np.log10(np.loadtxt(fileobs,skiprows=5, usecols = (2), unpack = True))
            yhist_err = np.abs(np.log10(np.loadtxt(fileobs,skiprows=5, usecols = (3), unpack = True)))

            plt.plot(hist_obs, yhist_obs, 'o', label = 'moustakas_2013_z' + str(iz))
            #plt.errorbar(hist_obs, yhist_obs, yhist_err,fmt='none', elinewidth=None, ecolor = 'red')

            # Observations Ilbert
            fileobs = '../Obs_Data/smf/muzzin_ilbert_z0.2_z0.5.smf'
            mass_low = np.loadtxt(fileobs, skiprows=4, usecols=(0), unpack=True) - np.log10(h)
            mass_high = np.loadtxt(fileobs, skiprows=4, usecols=(1), unpack=True) - np.log10(h)

            dm_obs = mass_high[1]-mass_low[1] # 0.25
            #print(dm_obs) # = 0.25
            hist_obs = mass_high - 0.5 * dm_obs
            yhist_obs = np.loadtxt(fileobs,skiprows=4, usecols = (2), unpack = True)
            yhist_err = np.loadtxt(fileobs,skiprows=4, usecols = (3), unpack = True)

            plt.plot(hist_obs, yhist_obs, 'o', label = 'ilbert_2013_z' + str(iz))
            #plt.errorbar(hist_obs, yhist_obs, yhist_err,fmt='none', elinewidth=None, ecolor = 'purple')


        if iz==1:
            # Observations Moustakas
            fileobs = '../Obs_Data/smf/moustakas_z0.80_z1.00.smf'
            mass_low = np.loadtxt(fileobs, skiprows=4, usecols=(0), unpack=True) - np.log10(h)
            mass_high = np.loadtxt(fileobs, skiprows=4, usecols=(1), unpack=True) - np.log10(h)

            dm_obs = mass_high[1]-mass_low[1] # 0.25
            #print(dm_obs) # = 0.25
            hist_obs = mass_high - 0.5 * dm_obs
            yhist_obs = np.loadtxt(fileobs,skiprows=4, usecols = (2), unpack = True)
            yhist_err = np.abs(np.log10(np.loadtxt(fileobs,skiprows=5, usecols = (3), unpack = True)))
            plt.plot(hist_obs, yhist_obs, 'o', label = 'moustakas_2013_z' + str(iz))
            #plt.errorbar(hist_obs, yhist_obs, yhist_err,fmt='none', elinewidth=None, ecolor = 'orange')

            # Observations Ilbert
            fileobs = '../Obs_Data/smf/muzzin_ilbert_z0.5_z1.1.smf'
            mass_low = np.loadtxt(fileobs, skiprows=4, usecols=(0), unpack=True) - np.log10(h)
            mass_high = np.loadtxt(fileobs, skiprows=4, usecols=(1), unpack=True) - np.log10(h)

            dm_obs = mass_high[1]-mass_low[1] # 0.25
            #print(dm_obs) # = 0.25
            hist_obs = mass_high - 0.5 * dm_obs
            yhist_obs = np.loadtxt(fileobs,skiprows=4, usecols = (2), unpack = True)
            yhist_err = np.loadtxt(fileobs,skiprows=4, usecols = (3), unpack = True)

            plt.plot(hist_obs, yhist_obs, 'o', label = 'ilbert_2013_z' + str(iz))
            #plt.errorbar(hist_obs, yhist_obs, yhist_err,fmt='none', elinewidth=None, ecolor = 'limegreen')


        if iz==2:
            # Observations Henriques
            fileobs='../Obs_Data/smf/henriques_2014_z2_cha.txt'
            mass_low = np.loadtxt(fileobs, skiprows=5, usecols=(0), unpack=True) - np.log10(h)
            mass_high = np.loadtxt(fileobs, skiprows=5, usecols=(1), unpack=True) - np.log10(h)

            dm_obs = mass_high[1]-mass_low[1] # 0.25
            #print(dm_obs) # = 0.25
            hist_obs = mass_high - 0.5 * dm_obs
            yhist_obs = np.log10(np.loadtxt(fileobs, skiprows=5, usecols=(2), unpack=True))#*dm)
            yhist_err = np.abs(np.log10(np.loadtxt(fileobs, skiprows=5, usecols=(3), unpack=True)))#*dm)

            plt.plot(hist_obs, yhist_obs, 'o', label = 'henriques_2014_z' + str(iz))
            #plt.errorbar(hist_obs, yhist_obs, yhist_err,fmt='none', elinewidth=None, ecolor = 'orange')


            # Observations Ilbert
            fileobs = '../Obs_Data/smf/muzzin_ilbert_z2.0_z2.5.smf'
            mass_low = np.loadtxt(fileobs, skiprows=4, usecols=(0), unpack=True) - np.log10(h)
            mass_high = np.loadtxt(fileobs, skiprows=4, usecols=(1), unpack=True) - np.log10(h)

            dm_obs = mass_high[1]-mass_low[1] # 0.25
            #print(dm_obs) # = 0.25
            hist_obs = mass_high - 0.5 * dm_obs
            yhist_obs = np.loadtxt(fileobs,skiprows=4, usecols = (2), unpack = True)
            yhist_err = np.loadtxt(fileobs,skiprows=4, usecols = (3), unpack = True)

            plt.plot(hist_obs, yhist_obs, 'o', label = 'ilbert_2013_z' + str(iz))
            #plt.errorbar(hist_obs, yhist_obs, yhist_err,fmt='none', elinewidth=None, ecolor = 'limegreen')

        plt.legend()
        plt.savefig('../plots/lm_z'+ str(iz)+'.pdf')
        plt.show()

plot_lm(['shark'],0.25, z=[0,1,2])
#exit()

def plot_sfr(model, dsfr,z=[0,1,2]):
    """Plot the sfrf for the different models and compare them all with the observations in the same plot for different redshifts.

    Args:
        model (array of str): name of the models
        dsfr (float): bin of sfr for the model's histograms
        z (array of float or int): the different redshifts. Defaults to [0,1,2].
    """
    h = 0.6774
    hobs = 0.71

    for iz in z:
        fig,ax=plt.subplots()
        for mod in model:
            if mod == 'shark':
                if iz == 0:
                    filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/128/multiple_batches/galaxies.hdf5'
                if iz == 1:
                    filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/97/multiple_batches/galaxies.hdf5'
                if iz == 2:
                    filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/78/multiple_batches/galaxies.hdf5'

            if mod == 'galform':
                filename = '../galform_output/galaxies.hdf5'

            with h5py.File(filename,"r") as f:
                if mod == 'shark':
                    volume = 1.56250E+07 #f['run_info/lbox'][()]**3 # [Mpc/h]³
                    #mres = f['run_info/particle_mass'][()] # [Msun/h] #! Es 0, mirar por qué
                    mres = 9.97e9 # [Msun/h]
                    mvir = f['galaxies/mvir_hosthalo'][:] # Dark matter mass of the host halo in which this galaxy resides [Msun/h]
                    sfr = f['galaxies/sfr_burst'][:] + f['galaxies/sfr_disk'][:] #[Msun/Gyr/h]
                    mstars = f['galaxies/mstars_bulge'][:] + f['galaxies/mstars_disk'][:] #[Msun/h]

                    #print(len(sfr))

                if mod == 'galform':
                    volume = 1.56250E+07 #1e9 # [Mpc/h]³ # more used-parameters
                    mres = 9.97e9 # [Msun/h]

                    sfr = f['Output001/mstardot'][:] + f['Output001/mstardot_burst'][:] # disk and bulge respectively [Msun/Gyr/h]
                    mstars = f['Output001/mstars_bulge'][:] + f['Output001/mstars_disk'][:]  #[Msun/h]

                    #print('works', len(sfr))

                ind = np.where(mvir>mres*20)
                lsfr = np.log10(sfr[ind]) - np.log10(h) - 9 #[Msun/yr]

                sfrmin = lsfr.min()
                sfrmax = lsfr.max()
                edges = np.array(np.arange(sfrmin, sfrmax+dsfr, dsfr)) #from sfrmin to sfrmax with a dex bin
                hist = edges[1:]-0.5*dsfr
                nsfr = np.zeros(shape=len(hist))

                iline=0

                for starfr in lsfr: # in log

                    for ie, edge in enumerate(edges[:-1]):
                        if starfr >= edge and starfr<edges[ie+1]:
                            nsfr[ie] = nsfr[ie] + 1

                    iline+=1

                #print(sum(nm))

                f.close

            yhist = np.log10(nsfr) - np.log10(dsfr) - np.log10(volume) + 3*np.log10(h)

            ax.plot(hist, yhist, label = mod)
            plt.title('Star Formation Rate Function z = ' + str(iz))
            plt.ylabel('log$_{10}$ ($\\Phi$ [Mpc$^{-3}$ dex$^{-1}$])')
            plt.xlabel('log$_{10}$(SFR[M$_{\\odot}$ yr$^{-1}$ ])')

        ### OBSERVATIONS ###
        # [h^-2Msun] If one wants to plot quantities without h in the units, then both observations and simulations should be divided by the h^2 of the simulation.
        if iz==0:
            fileobs='../Obs_Data/sfrf/gruppioni_2015_z0.0-0.3_cha.txt'

        if iz==1:
            fileobs='../Obs_Data/sfrf/gruppioni_2015_z0.8-1.0_cha.txt'

        if iz==2:
            fileobs='../Obs_Data/sfrf/gruppioni_2015_z2.0-2.5_cha.txt'
        sfr_low = np.loadtxt(fileobs, skiprows=3, usecols=(0), unpack=True) # log(SFR/(Msun/yr))
        sfr_high = np.loadtxt(fileobs, skiprows=3, usecols=(1), unpack=True) # log(SFR/(Msun/yr))

        dsfr_obs = sfr_high[1]-sfr_low[1]
        #print(dm_obs) # = 0.25
        hist_obs = sfr_high - 0.5 * dsfr_obs
        yhist_obs = np.loadtxt(fileobs, skiprows=3, usecols=(2), unpack=True)
        yhist_err = np.loadtxt(fileobs, skiprows=3, usecols=(3), unpack=True)

        plt.plot(hist_obs, yhist_obs, 'o', label = 'gruppioni_2015_z~' + str(iz))
        plt.errorbar(hist_obs, yhist_obs, yhist_err,fmt='none', elinewidth=None, ecolor = 'orange')

        plt.xlim(-2,5)
        plt.legend()
        plt.savefig('../plots/sfr_z'+ str(iz)+'.pdf')
        plt.show()

plot_sfr(['shark',], 0.25, z = [0,1,2])
exit()


def plot_BHSM(model, dbulge, dbh, z=0):
    """Plot the BHSM relation for the different models and compare them all with the observations in the same plot for different redshifts.

    Args:
        model (array of str): name of the models
        z (array of float or int): the different redshifts. Defaults to [0,1,2].
    """
    h = 0.6774

    fig,ax=plt.subplots()
    for mod in model:
        if mod == 'shark':
            filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/128/multiple_batches/galaxies.hdf5'
        if mod == 'galform':
            filename = '../galform_output/galaxies.hdf5'

        with h5py.File(filename,"r") as f:
            if mod == 'shark':
                volume = 1.56250E+07 #f['run_info/lbox'][()]**3 # [Mpc/h]³
                #mres = f['run_info/particle_mass'][()] # [Msun/h] #! Es 0, mirar por qué
                mres = 9.97e9 # [Msun/h]
                mvir = f['galaxies/mvir_hosthalo'][:] # Dark matter mass of the host halo in which this galaxy resides [Msun/h]

                mbh = f['galaxies/m_bh'][:] #[Msun/h]
                m_bulge = f['galaxies/mstars_bulge'][:]

                #print(len(sfr))

            if mod == 'galform':
                volume = 1.56250E+07 #1e9 # [Mpc/h]³ # more used-parameters
                mres = 9.97e9 # [Msun/h]
                print('no data yet')
                exit()

                #print('works', len(sfr))

            #mmin = np.log10(mres)
            #indmass = np.where(np.log10(mbh)>mmin)

            #lmbh = np.log10(mbh[indmass])
            #lm_bulge = np.log10(m_bulge[indmass])
            ind = np.where(mvir>mres*20)
            lm_bulge = np.log10(m_bulge[ind]) - np.log10(h) #log[Msun]
            lmbh = np.log10(mbh[ind]) - np.log10(h) #log[Msun]

            f.close
        ind = np.where(lm_bulge!=float('-inf'))
        #print(ind)
        #exit()
        lm_bulge = lm_bulge[ind]
        lmbh = lmbh[ind]
        bulge_min = lm_bulge.min()
        bulge_max = lm_bulge.max()
        print(bulge_min,bulge_max)

        bulge_edges = np.array(np.arange(bulge_min, bulge_max+dbulge, dbulge)) #from sfrmin to sfrmax with a dex bin
        
        bh_min = lmbh.min()
        bh_max = lmbh.max()
        bh_edges = np.array(np.arange(bh_min, bh_max+dbh, dbh)) #from sfrmin to sfrmax with a dex bin
        print(bh_min, bh_max)
        H, xedges, yedges = np.histogram2d(lm_bulge, lmbh, bins=([bulge_edges,bh_edges]))
        ind = np.where(H>0)
        median = np.median(H[ind])
        percentile_25 = np.percentile(H[ind], 25)
        percentile_75 = np.percentile(H[ind], 75)
        percentiles = np.percentile(H[ind], [25,50,75])
        print(percentiles)
        #exit()
        #ax.scatter(lm_bulge, lmbh, label = mod)

        #plt.contourf(xedges[:-1], yedges[:-1], H.T, levels=20, cmap='viridis')
        # Plot contours for median and percentiles
        #contour = plt.contour(xedges[:-1], yedges[:-1], H.T, levels=percentiles, colors='black')
        #plt.clabel(contour, fmt='%1.2f%%', inline=True)
        #X, Y = np.meshgrid(xedges, yedges)
        #ax.pcolormesh(X, Y, H)
        #plt.contourf(xedges[:-1], yedges[:-1], H.T, levels=20, cmap='viridis')
        plt.contour(xedges[:-1], yedges[:-1], H.T, levels=[median], colors='r', linestyles='solid', linewidths=2)
        plt.contour(xedges[:-1], yedges[:-1], H.T, levels=[percentile_25, percentile_75], colors='b', linestyles='dashed', linewidths=2)

        plt.title('BH-bulge mass relation z = 0')
        plt.ylabel('log$_{10}$ (M$_{BH}$/M$_{\\odot}$)')
        plt.xlabel('log$_{10}$ (M$_{bulge}$/M$_{\\odot}$)')

        ### OBSERVATIONS ###
        # McConnel
        fileobs='../Obs_Data/bhsm/McConnell_Ma_2013_ascii.txt'

        mbh_obs = np.log10(np.loadtxt(fileobs, skiprows=17, usecols=(2), unpack=True)) # [Msun]
        m_bulge_obs = np.log10(np.loadtxt(fileobs, skiprows=17, usecols=(11), unpack=True)) # [Msun]

        plt.plot(m_bulge_obs, mbh_obs, 'o', label = 'McConnel+2013')

        plt.ylim(5,11)
        plt.xlim(8,13)
        plt.legend()
        plt.savefig('../plots/BHSM_z0.pdf')
        plt.show()

#plot_BHSM(['shark'], 0.25, 0.25)
#exit()
def plot_sizes(model,dmbulge, drbulge, dmdisk, drdisk,dmtotal, z=[0]):
    """Plot the size-stellar mass relation for disks and bulges for the different models and compare them all with the observations in the same plot for different redshifts.

    Args:
        model (array of str): name of the models
        z (array of float or int): the different redshifts. Defaults to [0,1,2].
    """

    h = 0.6774
    hobs = 0.71


    for iz in z:
        fig_bulge,ax_bulge=plt.subplots()
        plt.title('Size-stellar mass relation for bulges z = ' + str(iz))
        plt.ylabel('log$_{10}$ (r$_{*, bulge}$/ckpc)')
        plt.xlabel('log$_{10}$ (M$_{*, bulge}$/M$_{\\odot}$)')
        plt.xlim(8,12)
        plt.ylim(-0.5,2)
        fig_disk,ax_disk=plt.subplots()
        plt.title('Size-stellar mass relation for disks z = ' + str(iz))
        plt.ylabel('log$_{10}$ (r$_{*, disk}$/ckpc)')
        plt.xlabel('log$_{10}$ (M$_{*, disk}$/M$_{\\odot}$)')
        plt.xlim(8,12)
        plt.ylim(-0.5,2)
        
        fig_total,ax_total=plt.subplots()
        plt.title('Size-stellar mass relation z = ' + str(iz))
        plt.ylabel('log$_{10}$ (r$_{*, disk}$/ckpc)')
        plt.xlabel('log$_{10}$ (M$_{*}$/M$_{\\odot}$)')
        plt.xlim(8,12)
        plt.ylim(-0.2,1.6)
        for mod in model:
            if mod == 'shark':
                if iz == 0:
                    filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/128/multiple_batches/galaxies.hdf5'
                if iz == 1:
                    filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/97/multiple_batches/galaxies.hdf5'
                if iz == 2:
                    filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/78/multiple_batches/galaxies.hdf5'

            if mod == 'galform':
                filename = '../galform_output/galaxies.hdf5'

            with h5py.File(filename,"r") as f:
                if mod == 'shark':
                    volume = 1.56250E+07 #f['run_info/lbox'][()]**3 # [Mpc/h]³
                    #mres = f['run_info/particle_mass'][()] # [Msun/h] #! Es 0, mirar por qué
                    mres = 9.97e9 # [Msun/h]
                    mvir = f['galaxies/mvir_hosthalo'][:] # Dark matter mass of the host halo in which this galaxy resides [Msun/h]

                    rstar_bulge = f['galaxies/rstar_bulge'][:] #[cMpc/h]
                    rstar_disk = f['galaxies/rstar_disk'][:] #[cMpc/h]
                    mstars_bulge = f['galaxies/mstars_bulge'][:] #[Msun/h]
                    mstars_disk = f['galaxies/mstars_disk'][:] #[Msun/h]

                    #print(len(sfr))

                if mod == 'galform':
                    volume = 1.56250E+07 #1e9 # [Mpc/h]³ # more used-parameters
                    mres = 9.97e9 # [Msun/h]
                    print('no data yet')
                    exit()

                    #print('works', len(sfr))


                ind = np.where(mvir>mres*20)

                lrstar_bulge = np.log10(rstar_bulge[ind]) + 3 - np.log10(h) #[ckpc]
                lrstar_disk = np.log10(rstar_disk[ind]) + 3 - np.log10(h) #[ckpc]
                lmstars_bulge = np.log10(mstars_bulge[ind]) - np.log10(h) #[Msun]
                lmstars_disk = np.log10(mstars_disk[ind]) - np.log10(h) #[Msun]
                f.close

            ## BULGE

            ind = np.where(lrstar_bulge!=float('-inf'))
            #print(ind)
            #exit()
            lrstar_bulge = lrstar_bulge[ind]
            lmstar_bulge = lmstars_bulge[ind]

            # Mass
            mbulge_min = lmstar_bulge.min()
            mbulge_max = lmstar_bulge.max()

            mbulge_edges = np.array(np.arange(mbulge_min, mbulge_max+dmbulge, dmbulge))
            #print(mbulge_min,mbulge_max, mbulge_edges)

            # Radius
            rbulge_min = lrstar_bulge.min()
            rbulge_max = lrstar_bulge.max()
            rbulge_edges = np.array(np.arange(rbulge_min, rbulge_max+drbulge, drbulge))
            #print(rbulge_min, rbulge_max, rbulge_edges)


            H, xedges, yedges = np.histogram2d(lmstar_bulge, lrstar_bulge, bins=([mbulge_edges,rbulge_edges]))

            ind = np.where(H>0)
            median = np.median(H[ind])
            percentile_25 = np.percentile(H[ind], 25)
            percentile_75 = np.percentile(H[ind], 75)
            percentiles = np.percentile(H[ind], [25,50,75])
            #print(percentiles)

            ax_bulge.contour(xedges[:-1], yedges[:-1], H.T, levels=[median], colors='r', linestyles='solid', linewidths=2)
            ax_bulge.contour(xedges[:-1], yedges[:-1], H.T, levels=[percentile_25, percentile_75], colors='b', linestyles='dashed', linewidths=2)


            ## DISK

            ind = np.where(lrstar_disk!=float('-inf'))
            #print(ind)
            #exit()
            lrstar_disk = lrstar_disk[ind]
            lmstar_disk = lmstars_disk[ind]

            # Mass
            mdisk_min = lmstar_disk.min()
            mdisk_max = lmstar_disk.max()

            mdisk_edges = np.array(np.arange(mdisk_min, mdisk_max+dmdisk, dmdisk))
            #print(mbulge_min,mbulge_max, mbulge_edges)

            # Radius
            rdisk_min = lrstar_disk.min()
            rdisk_max = lrstar_disk.max()
            rdisk_edges = np.array(np.arange(rdisk_min, rdisk_max+drdisk, drdisk))
            #print(rbulge_min, rbulge_max, rbulge_edges)


            H, xedges, yedges = np.histogram2d(lmstar_disk, lrstar_disk, bins=([mdisk_edges,rdisk_edges]))

            ind = np.where(H>0)
            median = np.median(H[ind])
            percentile_25 = np.percentile(H[ind], 25)
            percentile_75 = np.percentile(H[ind], 75)
            percentiles = np.percentile(H[ind], [25,50,75])
            #print(percentiles)

            ax_disk.contour(xedges[:-1], yedges[:-1], H.T, levels=[median], colors='r', linestyles='solid', linewidths=2)
            ax_disk.contour(xedges[:-1], yedges[:-1], H.T, levels=[percentile_25, percentile_75], colors='b', linestyles='dashed', linewidths=2)

            ## TOTAL

            ind = np.where(lrstar_disk!=float('-inf'))
            #print(ind)
            #exit()
            lrstar_disk = lrstar_disk[ind]
            lmstar = lmstars_disk[ind] + lmstars_bulge[ind]
            ind = np.where(lmstar!=float('-inf'))
            lmstar = lmstar[ind]
            lrstar_disk = lrstar_disk[ind]

            # Mass
            m_min = lmstar.min()
            m_max = lmstar.max()

            m_edges = np.array(np.arange(m_min, m_max+dmtotal, dmtotal))
            #print(mbulge_min,mbulge_max, mbulge_edges)

            # Radius
            rdisk_min = lrstar_disk.min()
            rdisk_max = lrstar_disk.max()
            rdisk_edges = np.array(np.arange(rdisk_min, rdisk_max+drdisk, drdisk))
            #print(rbulge_min, rbulge_max, rbulge_edges)


            H, xedges, yedges = np.histogram2d(lmstar, lrstar_disk, bins=([m_edges,rdisk_edges]))

            ind = np.where(H>0)
            median = np.median(H[ind])
            percentile_25 = np.percentile(H[ind], 25)
            percentile_75 = np.percentile(H[ind], 75)
            percentiles = np.percentile(H[ind], [25,50,75])
            #print(percentiles)

            ax_total.contour(xedges[:-1], yedges[:-1], H.T, levels=[median], colors='r', linestyles='solid', linewidths=2)
            ax_total.contour(xedges[:-1], yedges[:-1], H.T, levels=[percentile_25, percentile_75], colors='b', linestyles='dashed', linewidths=2)

            #ax1.scatter(lmstars_bulge, lrstar_bulge, label = mod)
            #ax2.scatter(lmstars_disk, lrstar_disk, label = mod)



        """
        ### OBSERVATIONS ###
        if iz==0:
            fileobs='data/gruppioni_2015_z0.0-0.3_cha.txt'

        if iz==1:
            fileobs='data/gruppioni_2015_z0.8-1.0_cha.txt'

        if iz==2:
            fileobs='data/gruppioni_2015_z2.0-2.5_cha.txt'
        sfr_low = np.loadtxt(fileobs, skiprows=3, usecols=(0), unpack=True) # log(SFR/(Msun/yr))
        sfr_high = np.loadtxt(fileobs, skiprows=3, usecols=(1), unpack=True) # log(SFR/(Msun/yr))

        plt.plot(hist_obs, yhist_obs, 'o', label = 'gruppioni_2015_z~' + str(iz))
        """

        fig_bulge.legend()
        fig_bulge.legend()
        fig_disk.savefig('..plots/sizes_bulge_z'+ str(iz)+'.pdf')
        fig_disk.savefig('../plots/sizes_disk_z'+ str(iz)+'.pdf')
        plt.show()

plot_sizes(['shark'], 0.25, 0.25, 0.25, 0.25, 0.35)
exit()

def plot_himf(model, dgas,z=[0]):
    """Plot the himf for the different models and compare them all with the observations in the same plot for different redshifts.

    Args:
        model (array of str): name of the models
        dgas (float): bin of gas mass for the model's histograms
        z (array of float or int): the different redshifts. Defaults to [0,1,2].
    """
    h = 0.6774
    for iz in z:
        fig,ax=plt.subplots()
        for mod in model:
            if mod == 'shark':
                if iz == 0:
                    filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/128/multiple_batches/galaxies.hdf5'
                if iz == 1:
                    filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/97/multiple_batches/galaxies.hdf5'
                if iz == 2:
                    filename = '../shark_output/UNIT-PNG/subf+subl+dhalos/78/multiple_batches/galaxies.hdf5'

            if mod == 'galform':
                filename = '../galform_output/galaxies.hdf5'

            with h5py.File(filename,"r") as f:
                if mod == 'shark':
                    volume = 1.56250E+07 #f['run_info/lbox'][()]**3 # [Mpc/h]³ #! es 1e9 el 1.56250E+07 es el volumen de 1 subvolumen, 1e9/64=1.56250E+07
                    #mres = f['run_info/particle_mass'][()] # [Msun/h] #! Es 0, mirar por qué
                    mres = 9.97e9 # [Msun/h]
                    mvir = f['galaxies/mvir_hosthalo'][:] # Dark matter mass of the host halo in which this galaxy resides [Msun/h]

                    mstars = f['galaxies/mstars_bulge'][:] + f['galaxies/mstars_disk'][:] #[Msun/h]
                    #print(len(mstars))
                    mgas = f['galaxies/matom_bulge'][:] + f['galaxies/matom_disk'][:] #[Msun/h]


                if mod == 'galform':
                    volume = 1.56250E+07 #1e09 # [Mpc/h]³ # more used-parameters
                    mres = 9.97e9 # [Msun/h]

                    mstars = f['Output001/mstars_bulge'][:] + f['Output001/mstars_disk'][:] #[Msun/h]
                    #print(len(mstars)) #! Hay muy pocas, porque tengo galform corrido en 1 solo subvolumen


                ind = np.where(mvir>mres*20)

                lgas = np.log10(mgas[ind])

                gasmin = lgas.min()
                gasmax = lgas.max()


                edges = np.array(np.arange(gasmin, gasmax+dgas, dgas)) #from mmin to mmax with a dex bin
                hist = edges[1:]-0.5*dgas
                ngas = np.zeros(shape=len(hist))

                iline=0

                for hi in lgas: # in log

                    for ie, edge in enumerate(edges[:-1]):
                        if hi >= edge and hi<edges[ie+1]:
                            ngas[ie] = ngas[ie] + 1

                    iline+=1

                #print(sum(nm))

                f.close

            yhist = np.log10(ngas) - np.log10(dgas) - np.log10(volume)

            ax.plot(hist, yhist, label = mod)
            plt.title('HI Mass Function z = ' + str(iz))

            plt.ylabel('log$_{10}$ (dn/dlog (M$_{*}$)/h$^{3}$ Mpc$^{-3}$)')
            plt.xlabel('log$_{10}$(M$_{*}$/h$^{-1}$ M$_{\\odot}$)')


        ### OBSERVATIONS ###
        # Not yet
        plt.legend()
        plt.savefig('../plots/himf_z'+ str(iz)+'.pdf')
        plt.show()

plot_himf(['shark'], 0.25)



