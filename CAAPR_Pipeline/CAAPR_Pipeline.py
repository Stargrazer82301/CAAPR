# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import time
import warnings
import numbers
import random
import shutil
import numpy as np
import scipy.ndimage
import multiprocessing as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.io.fits
import congrid
import ChrisFuncs
import CAAPR_IO
import CAAPR_Aperture
import CAAPR_Photom
import CAAPR_AstroMagic
sys.path.append( str( os.path.join( os.path.split( os.path.dirname(os.path.abspath(__file__)) )[0], 'CAAPR_AstroMagic', 'PTS') ) )
from pts.magic.catalog.importer import CatalogImporter





# The main pipeline; the cutout-production, aperture-fitting, and actual photometry parts of the CAAPR process are called in here, as sub-pipelines
def PipelineMain(source_dict, bands_dict, kwargs_dict):
    source_start = time.time()
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Processing target '+source_dict['name']+'.'



    # Convert cutout request strings to booleans, as necessary
    for band in bands_dict.keys():
        if bands_dict[band]['make_cutout']=='True':
            bands_dict[band]['make_cutout']=True
        if bands_dict[band]['make_cutout']=='False':
            bands_dict[band]['make_cutout']=False



        # Now check if cutouts are necessary; if so, produce them
        if bands_dict[band]['make_cutout']==True:
            raise ValueError('If you want to produce a cutout, please set the \'make_cutout\' field of the band table to be your desired cutout width, in arcsec.')
        if bands_dict[band]['make_cutout']>0:
            if 'band_dir_original' in bands_dict[band].keys():
                bands_dict[band]['band_dir'] = bands_dict[band]['band_dir_original']
            band_cutout_dir = CAAPR_IO.Cutout(source_dict, bands_dict[band], kwargs_dict['output_dir_path'], kwargs_dict['temp_dir_path'])

        # Otherwise, check if it is possible to trim padding of no-coverage from edge of map
        elif bands_dict[band]['make_cutout']==False:
            band_cutout_dir = CAAPR_IO.UnpaddingCutout(source_dict, bands_dict[band], kwargs_dict['output_dir_path'], kwargs_dict['temp_dir_path'])

        # Update current row of bands table to reflect the path of the freshly-made cutout
        if band_cutout_dir!=None:
            bands_dict[band]['band_dir_original'] = bands_dict[band]['band_dir']
            bands_dict[band]['band_dir'] = band_cutout_dir



    # Report to user that they're being a little bit odd
    if kwargs_dict['fit_apertures']==False and kwargs_dict['do_photom']==False:
        print '['+source_dict['name']+'] So you don\'t want aperture fitting, nor do you want actual photometry to happen? Erm, okay.'



    # Check if star-subtraction is requested for any band; if so, commence catalogue pre-fetching
    if kwargs_dict['starsub']==True:
        star_sub_check = False
        for band in bands_dict.keys():
            if bands_dict[band]['remove_stars']==True:
                star_sub_check = True
        if star_sub_check==True:
            if os.path.exists(os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic')):
                shutil.rmtree(os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic'))
            os.mkdir(os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic'))
            source_dict['pre_catalogue'] = CAAPR_AstroMagic.PreCatalogue(source_dict, bands_dict, kwargs_dict)



    # If aperture file not provided, commence aperture-fitting sub-pipeline
    if kwargs_dict['fit_apertures']==True:

        # In standard operation, process multiple sources in parallel
        aperture_start = time.time()
        aperture_output_list = []
        if kwargs_dict['parallel']==True:
            bands_dict_keys = bands_dict.keys()
            random.shuffle(bands_dict_keys)
            pool = mp.Pool(processes=kwargs_dict['n_proc'])
            for band in bands_dict_keys:
                aperture_output_list.append( pool.apply_async( CAAPR_Aperture.PipelineAperture, args=(source_dict, bands_dict[band], kwargs_dict) ) )
            pool.close()
            pool.join()
            del(pool)
            aperture_list = [output.get() for output in aperture_output_list if output.successful()==True]
            aperture_list = [aperture for aperture in aperture_list if aperture!=None]

        # If parallelisation is disabled, process sources one-at-a-time
        elif kwargs_dict['parallel']==False:
            for band in bands_dict.keys():
                aperture_output_list.append( CAAPR_Aperture.PipelineAperture(source_dict, bands_dict[band], kwargs_dict) )
                aperture_list = [output for output in aperture_output_list if output!=None]

        # Combine all fitted apertures to produce amalgam aperture
        aperture_combined = CAAPR_Aperture.CombineAperture(aperture_list, source_dict, kwargs_dict)

        # Record aperture properties to file
        aperture_string = str([ source_dict['name'], aperture_combined[0], aperture_combined[1], aperture_combined[2] ])#'name','semimaj_arcsec,axial_ratio,pos_angle\n'
        aperture_string = aperture_string.replace('[','').replace(']','').replace(' ','').replace('\'','')+'\n'
        aperture_table_file = open( kwargs_dict['aperture_table_path'], 'a')
        aperture_table_file.write(aperture_string)
        aperture_table_file.close()

        # Create grid of thumbnail images
        if kwargs_dict['thumbnails']==True:
            CAAPR_IO.ApertureThumbGrid(source_dict, bands_dict, kwargs_dict, aperture_list, aperture_combined)

        # Report time taken to fit apertures, and tidy up
        if kwargs_dict['verbose']: print '['+source_dict['name']+'] Time taken performing aperture fitting: '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(time.time()-aperture_start,4))+' seconds.'
        gc.collect()



    # Commence actual photometry sub-pipeline
    if kwargs_dict['do_photom']==True:

        # Handle problem where the user hasn't provided an aperture file, but also hasn't told CAAPR to fit its own apertures.
        if kwargs_dict['aperture_table_path']==False and kwargs_dict['fit_apertures']==False:
            raise ValueError('User has requested no aperture-fitting, and no photometry!')

        # In standard operation, process multiple sources in parallel
        photom_start = time.time()
        photom_output_list = []
        if kwargs_dict['parallel']==True:
            bands_dict_keys = bands_dict.keys()
            random.shuffle(bands_dict_keys)
            pool = mp.Pool(processes=kwargs_dict['n_proc'])
            for band in bands_dict_keys:
                photom_output_list.append( pool.apply_async( CAAPR_Photom.PipelinePhotom, args=(source_dict, bands_dict[band], kwargs_dict) ) )
            pool.close()
            pool.join()
            photom_list = [output.get() for output in photom_output_list if output.successful()==True]
            photom_list = [photom for photom in photom_list if photom!=None]

        # If parallelisation is disabled, process sources one-at-a-time
        elif kwargs_dict['parallel']==False:
            for band in bands_dict.keys():
                photom_output_list.append( CAAPR_Photom.PipelinePhotom(source_dict, bands_dict[band], kwargs_dict) )
                photom_list = [output for output in photom_output_list if output!=None]

        # Use band input table to establish order in which to put bands in results file, handling case of only one band
        photom_string = source_dict['name']
        bands_table = np.genfromtxt(kwargs_dict['bands_table_path'], delimiter=',', names=True, dtype=None)
        bands_list = bands_table['band_name']
        if bands_list.shape==():
            bands_list = [bands_list.tolist()]
        for band_name in bands_list:
            for photom_entry in photom_list:
                if photom_entry['band_name']==band_name:
                    photom_string += ','+str(photom_entry['ap_sum'])+','+str(photom_entry['ap_error'])
        photom_string += '\n'

        # Record photometry results to file
        photom_table_file = open( kwargs_dict['photom_table_path'], 'a')
        photom_table_file.write(photom_string)
        photom_table_file.close()

        # Create grid of thumbnail images
        if kwargs_dict['thumbnails']==True:
            CAAPR_IO.PhotomThumbGrid(source_dict, bands_dict, kwargs_dict)

        # Report time taken to do photometry, and tidy up
        if kwargs_dict['verbose']: print '['+source_dict['name']+'] Time taken performing actual photometry: '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(time.time()-photom_start,4))+' seconds.'



    # Tidy up temporary files
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Total time taken for souce: '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(time.time()-source_start,4))+' seconds.'
    if kwargs_dict['thumbnails']==True: [os.remove(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',processed_map)) for processed_map in os.listdir(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps')) if '.fits' in processed_map]
    if os.path.exists(os.path.join(kwargs_dict['temp_dir_path'],'Cutouts',source_dict['name'])):
        shutil.rmtree(os.path.join(kwargs_dict['temp_dir_path'],'Cutouts',source_dict['name']))
    if os.path.exists(os.path.join(kwargs_dict['temp_dir_path'],'AstroMagic')):
        shutil.rmtree(os.path.join(kwargs_dict['temp_dir_path'],'AstroMagic'))
    gc.collect()





# Define function that determines preliminary map values
def MapPrelim(pod, source_dict, band_dict, verbose=False):
    verbose = pod['verbose']



    # Check if x & y pixel sizes are meaningfully different. If so, panic; else, treat as same
    pix_size = 3600.0 * pod['in_wcs'].wcs.cdelt
    if float(abs(pix_size.max()))/float(abs(pix_size.min()))>(1+1E-3):
        raise ValueError('The x pixel size if noticably different from the y pixel size.')
    else:
        pod['pix_arcsec'] = float(np.mean(np.abs(pix_size)))

    # Determine source position in cutout in ij coordinates, and size of cutout
    centre_xy = pod['in_wcs'].wcs_world2pix( np.array([[ source_dict['ra'], source_dict['dec'] ]]), 0 )
    pod['centre_i'], pod['centre_j'] = float(centre_xy[0][1]), float(centre_xy[0][0])
    pod['box_rad'] = int( round( float(pod['cutout'].shape[0]) * 0.5 ) )

    # Determine beam size in pixels; if beam size not given, then assume map is Nyquist sampled (ie, 2.355 pixels ber beam)
    if isinstance(band_dict['beam_arcsec'], numbers.Number):
        pod['beam_pix'] = float(band_dict['beam_arcsec']) / pod['pix_arcsec']
    else:
        pod['beam_pix'] = pod['pix_arcsec'] * 2.355

    # Check if current source lies within bounds of map; if not, fai and return)
    if pod['centre_i']<0 or pod['centre_i']>(pod['cutout'].shape)[0] or pod['centre_j']<0 or pod['centre_j']>(pod['cutout'].shape)[1]:
        pod['within_bounds'] = False
        if verbose: print '['+pod['id']+'] Target not within bounds of map.'
    else:
        pod['within_bounds'] = True

    # Return pod
    return pod





# Define function that fits and subtracts polynomial background filter from map
def PolySub(pod, mask_semimaj_pix, mask_axial_ratio, mask_angle, poly_order=5, cutoff_sigma=2.0):
    if pod['verbose']: print '['+pod['id']+'] Determining if (and how) background is significnatly variable.'



    # Define Keflavich function to downsample an array
    def Downsample(myarr,factor,estimator=np.nanmean):
        ys,xs = myarr.shape
        crarr = myarr[:ys-(ys % int(factor)),:xs-(xs % int(factor))]
        dsarr = estimator( np.concatenate([[crarr[i::factor,j::factor]
            for i in range(factor)]
            for j in range(factor)]), axis=0)
        return dsarr



    # If image has pixels smaller than some limit, downsample image to improve processing time
    pix_size = pod['pix_arcsec']
    pix_size_limit = 2.0
    if pix_size<pix_size_limit:
        downsample_factor = int(np.ceil(pix_size_limit/pix_size))
    else:
        downsample_factor = 1
    image_ds = Downsample(pod['cutout'], downsample_factor)

    # Downsample related values accordingly
    mask_semimaj_pix = mask_semimaj_pix / downsample_factor
    centre_i = int(round(float((0.5*pod['centre_i'])-1.0)))
    centre_j = int(round(float((0.5*pod['centre_j'])-1.0)))



    # Find cutoff for excluding bright pixels by sigma-clipping map
    clip_value = ChrisFuncs.SigmaClip(image_ds, tolerance=0.01, sigma_thresh=3.0, median=True)
    noise_value = clip_value[0]
    field_value = clip_value[1]
    cutoff = field_value + ( cutoff_sigma * noise_value )



    # Mask all image pixels in masking region around source
    image_masked = image_ds.copy()
    ellipse_mask = ChrisFuncs.Photom.EllipseMask(image_ds, mask_semimaj_pix, mask_axial_ratio, mask_angle, centre_i, centre_j)
    image_masked[ np.where( ellipse_mask==1 ) ] = np.nan

    # Mask all image pixels identified as being high SNR
    image_masked[ np.where( image_masked>cutoff ) ] = np.nan

    # Use astropy to fit 2-dimensional polynomial to the image
    image_masked[ np.where( np.isnan(image_masked)==True ) ] = field_value
    poly_model = astropy.modeling.models.Polynomial2D(degree=poly_order)
    i_coords, j_coords = np.mgrid[:image_masked.shape[0], :image_masked.shape[1]]
    fitter = astropy.modeling.fitting.LevMarLSQFitter()
    i_coords = i_coords.flatten()
    j_coords = j_coords.flatten()
    image_flattened = image_masked.flatten()
    good = np.where(np.isnan(image_flattened)==False)
    i_coords = i_coords[good]
    j_coords = j_coords[good]
    image_flattened = image_flattened[good]
    fit = fitter(poly_model, i_coords, j_coords, image_flattened)

    # Create final polynomial filter (undoing downsampling using lorenzoriano GitHub script)
    i_coords, j_coords = np.mgrid[:image_ds.shape[0], :image_ds.shape[1]]
    poly_fit = fit(i_coords, j_coords)
    poly_full = scipy.ndimage.interpolation.zoom(poly_fit, [ float(pod['cutout'].shape[0])/float(poly_fit.shape[0]), float(pod['cutout'].shape[1])/float(poly_fit.shape[1]) ], mode='nearest') #poly_full = congrid.congrid(poly_fit, (pod['cutout'].shape[0], pod['cutout'].shape[1]), minusone=True)



    # Establish background variation before application of filter
    sigma_thresh = 3.0
    clip_in = ChrisFuncs.SigmaClip(pod['cutout'], tolerance=0.005, median=True, sigma_thresh=sigma_thresh)
    bg_in = pod['cutout'][ np.where( pod['cutout']<clip_in[1] ) ]
    spread_in = np.mean( np.abs( bg_in - clip_in[1] ) )

    # How much reduction in background variation there was due to application of the filter
    image_sub = pod['cutout'] - poly_full
    clip_sub = ChrisFuncs.SigmaClip(image_sub, tolerance=0.005, median=True, sigma_thresh=sigma_thresh)
    bg_sub = image_sub[ np.where( image_sub<clip_sub[1] ) ]
    spread_sub = np.mean( np.abs( bg_sub - clip_sub[1] ) )
    spread_diff = spread_in / spread_sub

    # If the filter made significant difference, apply to image and return it; otherwise, just return the unaltered map
    if spread_diff>1.1:
        if pod['verbose']: print '['+pod['id']+'] Background is significnatly variable; removing polynomial background fit.'
        pod['cutout_nopoly'] = pod['cutout'].copy()
        pod['cutout'] = image_sub
        pod['sky_poly'] = poly_model
    else:
        if pod['verbose']: print '['+pod['id']+'] Background is not significnatly variable; leaving image unaltered.'
        pod['sky_poly'] = False
    return pod





# Define function that predicts time until completion, and produces plot thereof
def TimeEst(time_list, total, output_dir_path, source_dict, kwargs_dict):
    time_list.append(time.time())
    time_est = ChrisFuncs.TimeEst(time_list, total, plot=True)
    time_remaining = time_est[0]
    time_file = open( os.path.join(output_dir_path,'Estimated_Completion_Time.txt'), 'w')
    time_file.write(time_remaining)
    time_file.close()
    time_fig = time_est[1]
    time_fig.savefig( os.path.join(output_dir_path,'Estimated_Completion_Time.png'), dpi=150  )
    time_fig.clf()
    plt.close('all')
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] CAAPR estimated completion at: '+time_remaining+'.'




