# Import smorgasbord
import sys
import os
import gc
import time
import random
#import warnings
#warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import multiprocessing as mp
import CAAPR
import CAAPR.CAAPR_IO
import CAAPR.CAAPR_Pipeline
import pdb




# Define the function that runs the CAAPR pipeline
def Run(bands_table_path = '../CAAPR_Example/CAAPR_Band_Table.csv',
        sources_table_path = '../CAAPR_Example/CAAPR_Source_Table.csv',
        output_dir_path = os.path.join(os.getcwd(),'CAAPR_Output'),
        temp_dir_path = os.path.join(os.getcwd(),'CAAPR_Temp'),
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
        save_images = False,
        debug = False,
        verbose = True,
        messy = False
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
                   'save_images': save_images,
                   'debug':debug,
                   'verbose':verbose,
                   'messy':messy}



    # Read in sources table, and convert into dictionary
    sources_dict = CAAPR.CAAPR_IO.SourcesDictFromCSV(sources_table_path)

    # Read in bands table, and convert into dictionary
    bands_dict = CAAPR.CAAPR_IO.BandsDictFromCSV(bands_table_path)

    # Prepare output directory
    CAAPR.CAAPR_IO.OutputDirPrepare(kwargs_dict)

    # Prepare temp directory, deleting any pre-existing directory at the specified location
    CAAPR.CAAPR_IO.TempDirPrepare(kwargs_dict)



    # Make inviolate copy of original band directories, to insure against over-writing when temp cutout directories are handled later
    for band in bands_dict.keys():
        bands_dict[band]['band_dir_inviolate'] = bands_dict[band]['band_dir']

    # Record timestamp
    kwargs_dict['timestamp'] = str(time.time()).replace('.','-')



    # If no aperture table file provided, and aperture-fitting is requested, create and prepare CSV file to store aperture dimensions for each source
    kwargs_dict = CAAPR.CAAPR_IO.ApertureTablePrepare(kwargs_dict)

    # If no photometry table path provided, and photometry is requested, create and prepare CSV file to store photometry output for each source
    kwargs_dict = CAAPR.CAAPR_IO.PhotomTablePrepare(kwargs_dict)




    # Randomise order of source dictionary keys (to "smooth out" average system resource usage)
    source_dict_keys = sources_dict.keys()
    random.shuffle(source_dict_keys)

    # Loop over each target source, processing in turn
    time_list = [time.time()]
    if verbose: print('[CAAPR] '+str(len(source_dict_keys))+' target objects to be processed.')
    for source in source_dict_keys:
        source_dict = sources_dict[source]
        CAAPR.CAAPR_Pipeline.PipelineMain(source_dict, bands_dict, kwargs_dict)

        # Estimate time until completions, and collect garbage
        CAAPR.CAAPR_Pipeline.TimeEst(time_list, len(source_dict_keys), output_dir_path, source_dict, kwargs_dict)
        gc.collect()





# Commence main task; generally you want to be calling CAAPR as a function, but it's useful to initiate a run this way for development and testing
if __name__ == "__main__":

    # Set parameters, and run function
    testing = True
    parallel = False
    starsub = True
    fit_apertures = True
    if fit_apertures==True:
        aperture_table_path = None
    elif fit_apertures==False:
        aperture_table_path = '../DustPedia/CAAPR_Aperture_Table_Test.csv'
    if testing:
        Run(temp_dir_path='/home/saruman/spx7cjc/DustPedia/CAAPR_Temp',
            n_proc=4,
            sources_table_path='../DustPedia/CAAPR_Source_Table_Test.csv',
            starsub=starsub,
            fit_apertures=fit_apertures,
            do_photom=False,
            aperture_table_path=aperture_table_path,
            parallel=parallel,
            debug=False,
            thumbnails=True)

        # Jubilate
        print('All done!')