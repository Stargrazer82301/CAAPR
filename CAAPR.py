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
import multiprocessing as mp

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
          photom_table_path = None,
          expansion_factor = 1.25,
          polysub = True,
          starsub = True,
          do_photom = True,
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
                   'parallel':parallel,
                   'n_proc':n_proc,
                   'thumbnails':thumbnails,
                   'verbose':verbose}

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
    if do_photom and thumbnails and not os.path.exists( os.path.join(output_dir_path,'Photometry_Thumbnails') ):
        os.mkdir( os.path.join(output_dir_path,'Photometry_Thumbnails') )

    # Prepare temp directory, deleting any pre-existing directory at the specified location
    if os.path.exists(temp_dir_path):
        shutil.rmtree(temp_dir_path)
    os.mkdir(temp_dir_path)
    if thumbnails==True:
        os.mkdir( os.path.join(temp_dir_path,'Processed_Maps') )
    os.mkdir(os.path.join(temp_dir_path, 'AstroMagic'))



    # Record timestamp
    timestamp = str(time.time()).replace('.','-')
    kwargs_dict['timestamp'] = timestamp

    # If no aperture table file provided, and aperture-fitting is requested, create and prepare CSV file to store aperture dimensions for each source
    if aperture_table_path==None and fit_apertures==True:
        aperture_table_path = os.path.join(output_dir_path,'CAAPR_Aperture_Table_'+timestamp+'.csv')
        kwargs_dict['aperture_table_path'] = aperture_table_path
        aperture_table_header = 'name,semimaj_arcsec,axial_ratio,pos_angle\n'
        aperture_table_file = open( aperture_table_path, 'a')
        aperture_table_file.write(aperture_table_header)
        aperture_table_file.close()

    # If no photometry table path provided, and photometry is requested, create and prepare CSV file to store photometry output for each source
    if photom_table_path==None and do_photom==True:
        photom_table_path = os.path.join(kwargs_dict['output_dir_path'],'CAAPR_Photom_Table_'+kwargs_dict['timestamp']+'.csv')
        kwargs_dict['photom_table_path'] = photom_table_path
        CAAPR_IO.PhotomTablePrepare(bands_dict, kwargs_dict)



    # Loop over each source, with completion time estimate
    source_dict_keys = sources_dict.keys()
    random.shuffle(source_dict_keys)
    time_list = [time.time()]
    for source in source_dict_keys:
        source_dict = sources_dict[source]
        CAAPR_Pipeline.PipelineMain(source_dict, bands_dict, kwargs_dict)
        time_list.append(time.time())
        time_remaining = ChrisFuncs.TimeEst(time_list, len(source_dict_keys))
        if verbose: print '['+source_dict['name']+'] CAAPR estimated completion at: '+time_remaining





# Commence main task; generally you want to be calling CAAPR as a function, but it's useful to initiate a run this way for development and testing
if __name__ == "__main__":

    # Run function
    testing = True
    parallel = True
    if testing:
        CAAPR(temp_dir_path='/home/saruman/spx7cjc/DustPedia/CAAPR_Temp', n_proc=10, polysub=True, starsub=False, sources_table_path='CAAPR_Source_Table.csv', fit_apertures=True, do_photom=False, aperture_table_path=None, parallel=parallel)

        # Jubilate
        print 'All done!'