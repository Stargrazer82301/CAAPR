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
import random
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp
import ChrisFuncs
import CAAPR_IO
import CAAPR_Pipeline
import CAAPR_Aperture
import CAAPR_Photom
import CAAPR_AstroMagic





def CAAPR(bands_table_path = 'CAAPR_Band_Table.csv',
          sources_table_path = 'CAAPR_Source_Table.csv',
          output_dir_path = 'CAAPR_Output',
          temp_dir_path = 'CAAPR_Temp',
          fit_apertures = True,
          aperture_table_path = None,#'CAAPR_Aperture_Table.csv',
          photom_table_path = None,
          expansion_factor = 1.25,
          polysub = True,
          starsub = True,
          do_photom = True,
          extinction_corr = True,
          parallel = True,
          n_proc = mp.cpu_count()-2,
          thumbnails = True,
          verbose = True
          ):



    # Create dictionary of kwarg values
    kwargs_dict = {'sources_table_path':sources_table_path,
                   'bands_table_path':bands_table_path,
                   'output_dir_path':output_dir_path,
                   'temp_dir_path':temp_dir_path,
                   'fit_apertures':fit_apertures,
                   'aperture_table_path':aperture_table_path,
                   'photom_table_path':photom_table_path,
                   'expansion_factor':expansion_factor,
                   'polysub':polysub,
                   'starsub':starsub,
                   'do_photom':do_photom,
                   'extinction_corr':extinction_corr,
                   'parallel':parallel,
                   'n_proc':n_proc,
                   'thumbnails':thumbnails,
                   'verbose':verbose}



    # Read in sources table and bands table, and convert into dictionaries
    sources_dict = CAAPR_IO.SourcesDictFromCSV(sources_table_path)
    bands_dict = CAAPR_IO.BandsDictFromCSV(bands_table_path)

    # Prepare output directory
    CAAPR_IO.OutputDirPrepare(kwargs_dict)

    # Prepare temp directory, deleting any pre-existing directory at the specified location
    CAAPR_IO.TempDirPrepare(kwargs_dict)



    # Make inviolate copy of original band directories, to insure against over-writing when temp cutout directories are handled later
    for band in bands_dict.keys():
        bands_dict[band]['band_dir_inviolate'] = bands_dict[band]['band_dir']

    # Record timestamp
    kwargs_dict['timestamp'] = str(time.time()).replace('.','-')



    # If no aperture table file provided, and aperture-fitting is requested, create and prepare CSV file to store aperture dimensions for each source
    kwargs_dict = CAAPR_IO.ApertureTablePrepare(kwargs_dict)

    # If no photometry table path provided, and photometry is requested, create and prepare CSV file to store photometry output for each source
    kwargs_dict = CAAPR_IO.PhotomTablePrepare(kwargs_dict)




    # Randomise order of source dictionary keys (to "smooth out" average system resource usage)
    source_dict_keys = sources_dict.keys()
    random.shuffle(source_dict_keys)

    # Loop over each target source, processing in turn
    time_list = [time.time()]
    if verbose: print '[CAAPR] '+str(len(source_dict_keys))+' target objects to be processed.'
    for source in source_dict_keys:
        source_dict = sources_dict[source]
        CAAPR_Pipeline.PipelineMain(source_dict, bands_dict, kwargs_dict)

        # Estimate time until completions, and collect garbage
        CAAPR_Pipeline.TimeEst(time_list, len(source_dict_keys), output_dir_path, source_dict, kwargs_dict)
        gc.collect()





# Commence main task; generally you want to be calling CAAPR as a function, but it's useful to initiate a run this way for development and testing
if __name__ == "__main__":

    # Run function
    testing = True
    parallel = False
    starsub = True
    aperture_table_path = None#'CAAPR_Aperture_Table_Test.csv'
    if testing:
        CAAPR(temp_dir_path='/home/saruman/spx7cjc/DustPedia/CAAPR_Temp', n_proc=7, sources_table_path='CAAPR_Source_Table_Test.csv', polysub=True, starsub=starsub, fit_apertures=True, do_photom=True, aperture_table_path=aperture_table_path, parallel=parallel)

        # Jubilate
        print 'All done!'