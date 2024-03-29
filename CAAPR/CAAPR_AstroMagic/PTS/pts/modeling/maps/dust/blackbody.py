#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.blackbody Contains the BlackBodyDustMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from ....evolve.engine import GAEngine, RawScoreCriteria
from ....evolve.genomes.list1d import G1DList
from ....evolve import mutators
from ....evolve import initializators
from ....evolve import constants

# Import the relevant PTS classes and modules
from ....core.tools.logging import log

# -----------------------------------------------------------------

k850 = 0.077

t1_min = 10.
t1_max = 30.
t1_step = 3.

t2_min = 30.
t2_max = 60.
t2_step = 10.

md_min = 4.
md_max = 6.
md_step = 0.02

ratio_min = 0.
ratio_max = 1.
ratio_guess = 0.5

# -----------------------------------------------------------------

class BlackBodyDustMapMaker(object):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(BlackBodyDustMapMaker, self).__init__()

        # -- Attributes --

    # -----------------------------------------------------------------

    def run(self, method="grid"):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # ...
        self.make_map(method)

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def make_map(self, method, plot=False):

        """
        This function ...
        :return:
        """

        data = np.genfromtxt("./eg_user_files/observations_stack.dat")
        dataID = np.genfromtxt("./eg_user_files/observations_stack.dat", usecols=0, dtype='string')

        wa = np.array([24, 100, 160, 250., 350., 500.])  # wavelengths in um

        # Loop over all pixels
        for i in range(len(data[:, 0])):

            ydata = [data[i, 30], data[i, 34], data[i, 36], data[i, 38], data[i, 40], data[i, 42]]
            yerr = [data[i, 31], data[i, 35], data[i, 37], data[i, 39], data[i, 41], data[i, 43]]
            D = data[i, 1] * 3 * 10 ** 5 / 67.30

            # Do the fit for this pixel
            if method == "grid": t1, t2, mdust, ratio = self.fit_grid(wa, ydata, yerr, D)
            elif method == "genetic": t1, t2, mdust, ratio = self.fit_genetic(wa, ydata, yerr, D)

            if plot:

                plt.plot(Mddist, Mdprob)
                plt.title(dataID[i])
                plt.xlabel('T_cold (K)')
                plt.show()

                print(Md_sav, T1_sav, T2_sav)

                plt.errorbar(wa, ydata, yerr=yerr, fmt='bo')
                x2 = np.arange(10., 600., 0.1)
                plt.loglog()
                plt.title(dataID[i])
                plt.ylim(0.001, 1)
                plt.plot(x2, two_blackbodies(x2, D, np.log10(Mdratio_sav) + Md_sav, T1_sav, np.log10(Mdratio_sav) + Md_sav, T2_sav), 'r', label="best fit")
                plt.xlabel('Wavelength (microns)')
                plt.ylabel('Flux (Jy)')
                plt.plot(x2, blackbody(x2, D, np.log10(Mdratio_sav) + Md_sav, T1_sav), ':', lw=2,
                         label="cold dust: logMd = %s, Tc= %s K " % (np.log10(Mdratio_sav) + Md_sav, T1_sav))
                plt.plot(x2, blackbody(x2, D, np.log10(Mdratio_sav) + Md_sav, T2_sav), ':', lw=2,
                         label="warm dust: logMd = %s, Tc= %s K " % (np.log10(Mdratio_sav) + Md_sav, T2_sav))
                plt.legend(loc=4)
                plt.show()

        plt.legend(frameon=False)
        plt.savefig('fit_bb.png')

    # -----------------------------------------------------------------

    def fit_grid(self, wa, ydata, yerr, D):

        """
        This function ...
        :return:
        """

        chi2 = float('inf')

        # Parameter ranges
        cold_temp_range = np.arange(t1_min, t1_max, t1_step)
        warm_temp_range = np.arange(t2_min, t2_max, t2_step)
        dust_mass_range = np.arange(md_min, md_max, md_step)

        Mdprob = np.zeros(len(dust_mass_range))

        # Best values
        cold_temp_best = None
        warm_temp_best = None
        dust_mass_best = None
        ratio_best = None

        ptot = 0
        # Loop over the temperatures and dust masses
        index = 0
        for dust_mass in dust_mass_range:
            for warm_temp in warm_temp_range:
                for cold_temp in cold_temp_range:

                    # Debugging
                    log.debug("Fitting dust ratio for a dust mass of " + str(dust_mass) + ", a warm temperature of " + str(warm_temp) + ", and a cold temperature of " + str(cold_temp) + " ...")

                    # Optimize
                    popt = minimize(leastsq, [ratio_guess], args=(dust_mass, wa, ydata, yerr, D, cold_temp, warm_temp), method='Nelder-Mead', options={'maxiter': 200})

                    Mdr_new = popt.x
                    chi2_new = leastsq(Mdr_new, dust_mass, wa, ydata, yerr, D, cold_temp, warm_temp)

                    if chi2 > chi2_new:

                        chi2 = chi2_new

                        # Set best parameters
                        ratio_best = Mdr_new
                        cold_temp_best = cold_temp
                        warm_temp_best = warm_temp
                        dust_mass_best = dust_mass

                        #print(Mddist[ii], Mdr_new, cold_temp_best, warm_temp_best)

                    prob = np.exp(-0.5 * chi2_new)
                    Mdprob[index] += prob
                    ptot += prob

                    index += 1

        Mdprob = Mdprob / ptot

        #print("tcold", Mddist[np.where(Mdprob == max(Mdprob))], percentiles(Mddist, Mdprob, 16), percentiles(Mddist, Mdprob, 50), percentiles(Mddist, Mdprob, 84))

        return cold_temp_best, warm_temp_best, dust_mass_best, ratio_best

    # -----------------------------------------------------------------

    def fit_genetic(self, wavelengths, ydata, yerr, D):

        """
        This function ...
        :param wavelengths:
        :param ydata:
        :param yerr:
        :param D:
        :return:
        """

        minima = [t1_min, t2_min, md_min, ratio_min]
        maxima = [t1_max, t2_max, md_max, ratio_max]

        # Create the first genome
        genome = G1DList(4)

        # Set genome options
        genome.setParams(minima=minima, maxima=maxima, bestrawscore=0.00, rounddecimal=2)
        genome.initializator.set(initializators.HeterogeneousListInitializerReal)
        # genome.mutator.set(mutators.HeterogeneousListMutatorRealRange)
        genome.mutator.set(mutators.HeterogeneousListMutatorRealGaussian)

        # Create the genetic algorithm engine
        engine = GAEngine(genome)

        # Set options for the engine
        engine.terminationCriteria.set(RawScoreCriteria)
        engine.setMinimax(constants.minimaxType["minimize"])
        engine.setGenerations(5)
        engine.setCrossoverRate(0.5)
        engine.setPopulationSize(100)
        engine.setMutationRate(0.5)

        # Initialize the genetic algorithm
        engine.initialize()

        ###

        engine.evolve()

        # Get best individual parameters
        best = engine.bestIndividual()
        best_t1 = best.genomeList[0]
        best_t2 = best.genomeList[1]
        best_md = best.genomeList[2]
        best_ratio = best.genomeList[3]

        # Return the best fitting parameters
        return best_t1, best_t2, best_md, best_ratio

# -----------------------------------------------------------------

def blackbody_base(lam, T):

    """
    Blackbody as a function of wavelength (um) and temperature (K).
    returns units of erg/s/cm^2/cm/Steradian
    """

    from scipy.constants import h, k, c
    lam = 1e-6 * lam  # convert to metres
    return 2. * h * c / (lam ** 3. * (np.exp(h * c / (lam * k * T)) - 1))

# -----------------------------------------------------------------

def blackbody(wa, D, Md, T1):

    """
    This function ...
    :param wa:
    :param D:
    :param Md:
    :param T1:
    :return:
    """

    kv = k850 * (850/wa)**2
    flux = kv * 10**Md/D**2 * blackbody_base(wa, T1)*2.08*10**11 #*1.98*10**30/9.5/10**32/10**12*10**26
    return flux

# -----------------------------------------------------------------

def two_blackbodies(wa, D, Md, T1, Md2, T2):

    """
    This function ...
    :param wa:
    :param D:
    :param Md:
    :param T1:
    :param Md2:
    :param T2:
    :return:
    """

    kv = k850 * (850/wa)**2
    flux = kv*10**Md/D**2*blackbody_base(wa, T1)*2.08*10**11+kv*10**Md2/D**2*blackbody_base(wa, T2)*2.08*10**11 #*1.98*10**30/9.5/10**32/10**12*10**26
    return flux

# -----------------------------------------------------------------

def leastsq(Dustratio, Md, wa, y, yerr, D, T1, T2):

    """
    This function ...
    :param Dustratio:
    :param Md:
    :param wa:
    :param y:
    :param yerr:
    :param D:
    :param T1:
    :param T2:
    :return:
    """

    som = 0
    y2 = two_blackbodies(wa, D, np.log10(Dustratio)+Md,T1,np.log10(1-(Dustratio))+Md,T2)

    for i in range(len(wa)):
        if i!=0 or y[i]<y2[i]:
            som += ((y[i]-y2[i])/yerr[i])**2

    return som

# -----------------------------------------------------------------

def percentiles(T1dist, T1prob, percentile):

    """
    This function ...
    :param T1dist:
    :param T1prob:
    :param percentile:
    :return:
    """

    percentile = percentile/100.
    perc = 0

    for ii in range(len(T1dist)):

        perc += T1prob[ii]
        if perc > percentile:
            return T1dist[ii]

# -----------------------------------------------------------------
