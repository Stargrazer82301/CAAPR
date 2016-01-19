# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import numpy as np
import multiprocessing as mp
#import matplotlib.pyplot as plt
#import astropy.io.fits
#import astropy.io.votable
#import astropy.wcs
#plt.ioff()

# Import ChrisFuncs and CAAPR submodules
import ChrisFuncs
import CAAPR_IO



def CAAPR(band_table_path='CAAPR_Band_Table.csv',
          source_table_path='CAAPR_Source_Table.csv',
          ):

    # Read in source tables
    band_table = np.genfromtxt(band_table_path, names=True, delimiter=',', dtype=None)
    source_table = np.genfromtxt(source_table_path, names=True, delimiter=',', dtype=None)

    # Commence multiprocesing
    output_list = []
    pool = mp.Pool()
    for i in range(0, len(source_table.shape)):
        source_dict = CAAPR_IO.MakeSourceDict(source_table, i)
        output_list.append( pool.apply_async( CAAPR.Pipeline, args=(source_dict,) ) )
    pool.close()
    pool.join()

    # Extract results
    output_list = [output.get() for output in output_list]



# Commence main task; generally you want to be calling CAAPR as a function, but it's useful to initiate a run this way for development and testing
if __name__ == "__main__":

    # Run function
    testing = True
    if testing:
        CAAPR()

        # Jubilate
        print 'All done!'