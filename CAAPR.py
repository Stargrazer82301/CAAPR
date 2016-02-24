# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import shutil
import time
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
          fit_apertures = True,
          aperture_table_path = None,
          parallel = True,
          n_cores = mp.cpu_count()-2,
          thumbnails = True,
          verbose = True
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
    if thumbnails==True:
        os.mkdir( os.path.join(temp_dir_path,'Thumbnail_Maps') )

    # Prepare CSV file to store aperture dimensions for each source
    timestamp = str(time.time()).replace('.','')
    if aperture_table_path==None:
        aperture_table_path = os.path.join(output_dir_path,'CAAPR_Aperture_Table_'+timestamp+'.csv')
    aperture_table_header = 'name,semimaj_arcsec,axial_ratio,pos_angle\n'
    aperture_table_file = open( aperture_table_path, 'a')
    aperture_table_file.write(aperture_table_header)
    aperture_table_file.close()

    # Loop over each source
    for source in sources_dict.keys():
        source_dict = sources_dict[source]
        CAAPR_Pipeline.PipelineMain(source_dict, bands_dict, output_dir_path, temp_dir_path, fit_apertures, aperture_table_path, parallel, n_cores, thumbnails, verbose)



# Commence main task; generally you want to be calling CAAPR as a function, but it's useful to initiate a run this way for development and testing
if __name__ == "__main__":

    # Run function
    testing = True
    parallel = True
    if testing:
        CAAPR(parallel=parallel)

        # Jubilate
        print 'All done!'