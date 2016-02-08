# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import shutil
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
import CAAPR_Pipeline



def CAAPR(bands_table_path = 'CAAPR_Band_Table.csv',
          sources_table_path = 'CAAPR_Source_Table.csv',
          output_dir_path = 'CAAPR_Output',
          temp_dir_path = 'CAAPR_Temp',
          n_cores = mp.cpu_count()-2,
          ):

    # Read in source table and band table CSVs, and convert into dictionaries
    sources_dict = CAAPR_IO.SourcesDictFromCSV(sources_table_path)
    bands_dict = CAAPR_IO.BandsDictFromCSV(bands_table_path)

    # Prepare output directory
    if os.path.exists(output_dir_path):
        print 'Warning: Output directory already exists; some files may be overridden'
    else:
        os.mkdir(output_dir_path)

    # Prepare temp directory, deleting any pre-existing directory at the specified location
    if os.path.exists(temp_dir_path):
        shutil.rmtree(temp_dir_path)
    os.mkdir(temp_dir_path)

    # Loop over each source
    for source in sources_dict.keys():
        source_dict = sources_dict[source]
        CAAPR_Pipeline.PipelineMain(source_dict, bands_dict, output_dir_path, temp_dir_path, n_cores)



# Commence main task; generally you want to be calling CAAPR as a function, but it's useful to initiate a run this way for development and testing
if __name__ == "__main__":

    # Run function
    testing = True
    if testing:
        CAAPR()

        # Jubilate
        print 'All done!'