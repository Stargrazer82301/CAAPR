# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import shutil
import psutil
import time
import resource
import warnings
warnings.filterwarnings('ignore')
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
          aperture_table_path = None,#'CAAPR_Aperture_Table.csv',
          expansion_factor = 1.25,
          do_photom = False,
          parallel = True,
          n_cores = mp.cpu_count()-4,
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
    if fit_apertures and thumbnails and not os.path.exists( os.path.join(output_dir_path,'Aperture_Fitting_Thumbnails') ):
        os.mkdir( os.path.join(output_dir_path,'Aperture_Fitting_Thumbnails') )

    # Prepare temp directory, deleting any pre-existing directory at the specified location
    if os.path.exists(temp_dir_path):
        shutil.rmtree(temp_dir_path)
    os.mkdir(temp_dir_path)
    if thumbnails==True:
        os.mkdir( os.path.join(temp_dir_path,'Processed_Maps') )
    os.mkdir(os.path.join(temp_dir_path, 'AstroMagic'))

    # If no aperture table file provided, create and prepare CSV file to store aperture dimensions for each source
    timestamp = str(time.time()).replace('.','')
    if aperture_table_path==None:
        aperture_table_path = os.path.join(output_dir_path,'CAAPR_Aperture_Table_'+timestamp+'.csv')
        aperture_table_header = 'name,semimaj_arcsec,axial_ratio,pos_angle\n'
        aperture_table_file = open( aperture_table_path, 'a')
        aperture_table_file.write(aperture_table_header)
        aperture_table_file.close()

    # Create dictionary of kwarg values
    kwargs_dict = {'sources_table_path':sources_table_path,
                   'bands_table_path':bands_table_path,
                   'output_dir_path':output_dir_path,
                   'temp_dir_path':temp_dir_path,
                   'fit_apertures':fit_apertures,
                   'aperture_table_path':aperture_table_path,
                   'expansion_factor':expansion_factor,
                   'do_photom':do_photom,
                   'parallel':parallel,
                   'n_cores':n_cores,
                   'thumbnails':thumbnails,
                   'verbose':verbose}

    # Loop over each source
    for source in sources_dict.keys():
        source_dict = sources_dict[source]
        CAAPR_Pipeline.PipelineMain(source_dict, bands_dict, kwargs_dict)#sources_table_path, bands_table_path, output_dir_path, temp_dir_path, fit_apertures, aperture_table_path, parallel, n_cores, thumbnails, verbose)



# Commence main task; generally you want to be calling CAAPR as a function, but it's useful to initiate a run this way for development and testing
if __name__ == "__main__":

    # Run function
    testing = True
    parallel = False
    if testing:
        CAAPR(temp_dir_path='/home/saruman/spx7cjc/DustPedia/CAAPR_Temp', parallel=parallel)

        # Jubilate
        print 'All done!'





"""
# Put in place RAM limit
if ram_limit!=False:
    resource.setrlimit(resource.RLIMIT_AS, ( int(float(ram_limit)*float(psutil.virtual_memory()[0])), int(float(ram_limit)*float(psutil.virtual_memory()[0])) ) )
"""