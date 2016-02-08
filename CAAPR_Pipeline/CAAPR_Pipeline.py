# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import numpy as np
import multiprocessing as mp
import CAAPR_IO
import CAAPR_Aperture



# The main pipeline; the cutout-production, aperture-fitting, and actual photometry parts of the CAAPR process are called in here, as sub-pipelines
def PipelineMain(source_dict, bands_dict, output_dir_path, temp_dir_path, n_cores):

    # Check if cutouts are necessary; if so, produce them
    for band in bands_dict.keys():
        if bands_dict[band]['make_cutout']==True:
            raise ValueError('If you want to produce a cutout, please set the \'make_cutout\' field of the band table to be your desired cutout width, in arcsec.')
        if bands_dict[band]['make_cutout']>0:
            band_cutout_path = CAAPR_IO.Cutout(source_dict, bands_dict[band], output_dir_path, temp_dir_path)

            # Update current row of bands table to reflect the path of the freshly-made cutout
            bands_dict[band]['band_path'] = band_cutout_path
            pdb.set_trace()

    # Commence aperture-fitting sub-pipeline, processing multiple sources in parallel
    output_list = []
    pool = mp.Pool()
    for band in bands_dict.keys():
        output_list.append( pool.apply_async( CAAPR_Aperture.PipelineAperture, args=(source_dict, bands_dict[band], output_dir_path, temp_dir_path,) ) )
    pool.close()
    pool.join()

    # Extract and sort results of aperture-fitting sub-pipeline
    output_list = [output.get() for output in output_list]


