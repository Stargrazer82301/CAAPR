# Import smorgasbord
import os
import sys
sys.path.append( str( os.path.join( os.path.split( os.path.dirname(os.path.abspath(__file__)) )[0], 'CAAPR', 'CAAPR_AstroMagic', 'PTS') ) )
import gc
import pdb
import time
import re
import copy
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
import ChrisFuncs
import ChrisFuncs.Photom
import ChrisFuncs.FromGitHub
import CAAPR





# The main pipeline; the cutout-production, aperture-fitting, and actual photometry parts of the CAAPR process are called in here, as sub-pipelines
def PipelineMain(source_dict, bands_dict, kwargs_dict):



    # Start timer, and check that the user has actually asked CAAPR to do something; if they haven't asked CAAPR to do anything at all, tell them that they're being a bit odd!
    source_start = time.time()
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Processing target '+source_dict['name']+'.'
    if kwargs_dict['fit_apertures']==False and kwargs_dict['do_photom']==False:
        print '['+source_dict['name']+'] So you don\'t want aperture fitting, nor do you want actual photometry to happen? Erm, okay.'



    # Check if any data actually exists for this source
    if SourcePrelim(source_dict, bands_dict, kwargs_dict)==False:
        return



    # Loop over bands for initial processing
    for band in bands_dict.keys():

        # Do basic initial handling of band parameters
        bands_dict[band] = BandInitiate(bands_dict[band])

        # Functiont hat checks if user has requested a cutout; and, if so, produces it
        bands_dict[band] = CAAPR.CAAPR_IO.Cutout(source_dict, bands_dict[band], kwargs_dict)

        # Function that check sif it is possible to trim padding of no-coverage from edge of map (if user hasn't specificed a particular cutout be made)
        bands_dict[band] = CAAPR.CAAPR_IO.UnpaddingCutout(source_dict, bands_dict[band], kwargs_dict)



    # Check if star-subtraction is requested for any band; if so, commence catalogue pre-fetching
    CAAPR.CAAPR_AstroMagic.PreCatalogue(source_dict, bands_dict, kwargs_dict)



    # If aperture file not provided, commence aperture-fitting sub-pipeline
    if kwargs_dict['fit_apertures']==True:

        # Process sources inside while loop, to catch 'missed' bands
        aperture_attempts = 0
        while aperture_attempts!='Success':

            # In standard operation, process multiple sources in parallel
            aperture_start = time.time()
            aperture_output_list = []
            if kwargs_dict['parallel']==True:
                bands_dict_keys = bands_dict.keys()
                random.shuffle(bands_dict_keys)
                pool = mp.Pool(processes=kwargs_dict['n_proc'])
                for band in bands_dict_keys:
                    aperture_output_list.append( pool.apply_async( CAAPR.CAAPR_Aperture.SubpipelineAperture, args=(source_dict, bands_dict[band], kwargs_dict) ) )
                pool.close()
                pool.join()
                del(pool)
                aperture_list = [output.get() for output in aperture_output_list if output.successful()==True]
                aperture_list = [aperture for aperture in aperture_list if aperture!=None]

            # If parallelisation is disabled, process sources one-at-a-time
            elif kwargs_dict['parallel']==False:
                for band in bands_dict.keys():
                    aperture_output_list.append( CAAPR.CAAPR_Aperture.SubpipelineAperture(source_dict, bands_dict[band], kwargs_dict) )
                    aperture_list = [output for output in aperture_output_list if output!=None]

            # Check that all photometry completed
            aperture_attempts = CAAPR.CAAPR_Aperture.ApertureCheck(aperture_attempts, aperture_output_list, source_dict, bands_dict, kwargs_dict)

        # Combine all fitted apertures to produce amalgam aperture
        aperture_combined = CAAPR.CAAPR_Aperture.CombineAperture(aperture_list, source_dict, kwargs_dict)

        # Record aperture properties to file
        CAAPR.CAAPR_IO.RecordAperture(aperture_combined, source_dict, kwargs_dict)

        # Prepare thumbnail images for bands excluded from aperture fitting
        CAAPR.CAAPR_Aperture.ExcludedThumb(source_dict, bands_dict, kwargs_dict, aperture_list, aperture_combined)

        # Create grid of thumbnail images
        CAAPR.CAAPR_IO.ApertureThumbGrid(source_dict, bands_dict, kwargs_dict, aperture_list, aperture_combined)

        # Report time taken to fit apertures, and tidy up garbage
        gc.collect()
        if kwargs_dict['verbose']: print '['+source_dict['name']+'] Time taken performing aperture fitting: '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(time.time()-aperture_start,4))+' seconds.'




    # Commence actual photometry sub-pipeline
    if kwargs_dict['do_photom']==True:

        # Handle problem where the user hasn't provided an aperture file, but also hasn't told CAAPR to fit its own apertures.
        if kwargs_dict['aperture_table_path']==False and kwargs_dict['fit_apertures']==False:
            raise Exception('User has requested no aperture-fitting, and no photometry!')

        # Process sources inside while loop, to catch 'missed' bands
        photom_attempts = 0
        while photom_attempts!='Complete':

                # In standard operation, process multiple sources in parallel
                photom_start = time.time()
                photom_output_list = []
                if kwargs_dict['parallel']==True:
                    bands_dict_keys = bands_dict.keys()
                    random.shuffle(bands_dict_keys)
                    pool = mp.Pool(processes=kwargs_dict['n_proc'])
                    for band in bands_dict_keys:
                        photom_output_list.append( pool.apply_async( CAAPR.CAAPR_Photom.SubpipelinePhotom, args=(source_dict, bands_dict[band], kwargs_dict) ) )
                    pool.close()
                    pool.join()
                    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Gathering parallel threads.'
                    photom_output_list = [output.get() for output in photom_output_list if output.successful()==True]
                    photom_list = [photom for photom in photom_output_list if photom!=None]

                # If parallelisation is disabled, process sources one-at-a-time
                elif kwargs_dict['parallel']==False:
                    for band in bands_dict.keys():
                        photom_output_list.append( CAAPR.CAAPR_Photom.SubpipelinePhotom(source_dict, bands_dict[band], kwargs_dict) )
                        photom_list = [photom for photom in photom_output_list if photom!=None]

                # Check that all photometry completed
                photom_attempts, photom_output_list = CAAPR.CAAPR_Photom.PhotomCheck(photom_attempts, photom_output_list, source_dict, bands_dict, kwargs_dict)



        # Record photometry results to file
        CAAPR.CAAPR_IO.RecordPhotom(photom_list, source_dict, bands_dict, kwargs_dict)

        # Prepare thumbnail images for bands excluded from photometry
        CAAPR.CAAPR_Photom.ExcludedThumb(source_dict, bands_dict, kwargs_dict)

        # Create grid of thumbnail images
        CAAPR.CAAPR_IO.PhotomThumbGrid(source_dict, bands_dict, kwargs_dict)

        # Report time taken to do photometry, and tidy up
        if kwargs_dict['verbose']: print '['+source_dict['name']+'] Time taken performing actual photometry: '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(time.time()-photom_start,4))+' seconds.'



    # Tidy up temporary files and paths
    bands_dict = PathTidy(source_dict, bands_dict, kwargs_dict)

    # Report time taken for source, and tidy up garbage
    gc.collect()
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Total time taken for souce: '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(time.time()-source_start,4))+' seconds.'
    if kwargs_dict['thumbnails']==True: [os.remove(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',processed_map)) for processed_map in os.listdir(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps')) if '.fits' in processed_map]





# Define function to check if data actually exists for any band for this source
def SourcePrelim(source_dict, bands_dict, kwargs_dict):



    # Check that any of the bands actually have data for this source
    kwargs_dict_copy = copy.deepcopy(kwargs_dict)
    kwargs_dict_copy['verbose'] = False    
    bands_check = []
    for band in bands_dict.keys():
        source_id = source_dict['name']+'_'+bands_dict[band]['band_name']
        in_fitspath, file_found = CAAPR.CAAPR_Pipeline.FilePrelim(source_dict, bands_dict[band], kwargs_dict_copy)
        bands_check.append(file_found)

    # Report to user if no data found
    if True not in bands_check:
        print '['+source_id+'] No data found in target directory for current source.'

        # Make null entries in tables, as necessary
        if kwargs_dict['fit_apertures']==True:
            null_aperture_combined = [np.NaN, np.NaN, np.NaN, np.NaN]
            CAAPR.CAAPR_IO.RecordAperture(null_aperture_combined, source_dict, kwargs_dict)
        if kwargs_dict['do_photom']==True:
            CAAPR.CAAPR_IO.RecordPhotom([], source_dict, bands_dict, kwargs_dict)

    # Return result
    if True not in bands_check:
        return False
    elif True in bands_check:
        return True





# Define function that does basic initial handling of band parameters
def BandInitiate(band_dict):



    # Make sure band has content
    if band_dict==None:
        return band_dict

    # Parse band cutout request, converting string to boolean if necessary
    if band_dict['make_cutout']=='True':
        band_dict['make_cutout']=True
    elif band_dict['make_cutout']=='False':
        band_dict['make_cutout']=False
    else:
        try:
            band_dict['make_cutout'] = float(band_dict['make_cutout'])
        except:
            raise Exception('Cutout request not understood; should either be False, or width of cutout in arcseconds.')

    # Reset band directory to inviolate value, to purge any holdovers from previous source
    band_dict['band_dir'] = band_dict['band_dir_inviolate']

    # Return band dict
    return band_dict




# Define function that performs preimilary checks of file type and location
def FilePrelim(source_dict, band_dict, kwargs_dict):



    # Determine whether the user is specificing a directroy full of FITS files in this band (in which case use standardised filename format), or just a single FITS file
    try:
        if os.path.isdir(band_dict['band_dir']):
            in_fitspath = os.path.join( band_dict['band_dir'], source_dict['name']+'_'+band_dict['band_name'] )
        elif os.path.isfile(band_dict['band_dir']):
            in_fitspath = os.path.join( band_dict['band_dir'] )
    except:
        pdb.set_trace()
        
    # Work out whether the file extension for FITS file in question is .fits or .fits.gz
    file_found = False
    try:
        if os.path.exists(in_fitspath+'.fits'):
            in_fitspath = in_fitspath+'.fits'
            file_found = True
        elif os.path.exists(in_fitspath+'.fits.gz'):
            in_fitspath = in_fitspath+'.fits.gz'
            file_found = True
    except:
        raise Exception('Path provided for band '+str(band_dict['band_name'])+' refers to neither a file nor a folder.')
    
    # Return file values
    return in_fitspath, file_found





# Initiate the pod (Photometry Organisation Dictionary)
def PodInitiate(in_fitspath, source_dict, band_dict, kwargs_dict):
    source_id = source_dict['name']+'_'+band_dict['band_name']
    if kwargs_dict['verbose']: print '['+source_id+'] Reading in FITS data.'



    # Read in FITS file in question
    in_fitsdata = astropy.io.fits.open(in_fitspath)
    in_image = in_fitsdata[0].data
    in_header = in_fitsdata[0].header
    in_fitsdata.close()
    in_wcs = astropy.wcs.WCS(in_header)
    in_fitspath_size = float(os.stat(in_fitspath).st_size)

    # Create the pod (Photometry Organisation Dictionary), which will bundle all the photometry data for this source & band into one dictionary to be passed between functions
    pod = {'in_fitspath':in_fitspath,
           'in_image':in_image,
           'in_header':in_header,
           'in_wcs':in_wcs,
           'cutout':in_image.copy(),
           'output_dir_path':kwargs_dict['output_dir_path'],
           'temp_dir_path':kwargs_dict['temp_dir_path'],
           'in_fitspath_size':in_fitspath_size,
           'id':source_id,
           'verbose':kwargs_dict['verbose']}

    # Return pod
    return pod





# Define function that determines preliminary map values
def MapPrelim(pod, source_dict, band_dict, verbose=False):
    if pod['verbose']: print '['+pod['id']+'] Determining properties of map.'



    # Check if x & y pixel sizes are meaningfully different. If so, panic; else, treat as same
    pix_size = 3600.0 * pod['in_wcs'].wcs.cdelt
    if float(abs(pix_size.max()))/float(abs(pix_size.min()))>(1+1E-3):
        raise Exception('The x pixel size if noticably different from the y pixel size.')
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
        if 'band_dir_inviolate' in band_dict.keys():
            band_dict['band_dir'] = band_dict['band_dir_inviolate']
        if pod['verbose']: print '['+pod['id']+'] Target not within bounds of map.'
    else:
        pod['within_bounds'] = True

    # Return pod
    return pod





# Define function that fits and subtracts polynomial background filter from map
def PolySub(pod, mask_semimaj_pix, mask_axial_ratio, mask_angle, poly_order=5, cutoff_sigma=2.0, instant_quit=False):
    if pod['verbose']: print '['+pod['id']+'] Determining if (and how) background is significnatly variable.'



    # If polynomial background subraction not wanted, immediately return everything unchanged
    if instant_quit:
         pod['sky_poly'] = False
         return pod



    # If image has pixels smaller than some limit, downsample image to improve processing time
    pix_size = pod['pix_arcsec']
    pix_size_limit = 2.0
    if pix_size<pix_size_limit:
        downsample_factor = int(np.ceil(pix_size_limit/pix_size))
    else:
        downsample_factor = 1
    image_ds = ChrisFuncs.Downsample(pod['cutout'], downsample_factor)

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

    # Use astropy to set up 2-dimensional polynomial to the image
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

    # Attempt polynomial fit; if insufficient data then skip onwards
    image_flattened = image_flattened[good]
    try:
        fit = fitter(poly_model, i_coords, j_coords, image_flattened)
    except:
        if pod['verbose']: print '['+pod['id']+'] Background is not significnatly variable; leaving image unaltered.'
        pod['sky_poly'] = False
        return pod

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




# Define function that tidies up folders and paths after completed processing a source
def PathTidy(source_dict, bands_dict, kwargs_dict):

    if os.path.exists(os.path.join(kwargs_dict['temp_dir_path'],'Cutouts',source_dict['name'])):
        shutil.rmtree(os.path.join(kwargs_dict['temp_dir_path'],'Cutouts',source_dict['name']))
    if os.path.exists(os.path.join(kwargs_dict['temp_dir_path'],'AstroMagic')):
        shutil.rmtree(os.path.join(kwargs_dict['temp_dir_path'],'AstroMagic'))

    # Set band directories to standard, not whatever temporary cutout directories may have been used for this source
    for band in bands_dict.keys():
        if bands_dict[band]==None:
            continue
        bands_dict[band]['band_dir'] = bands_dict[band]['band_dir_inviolate']



# Define function that predicts time until completion, and produces plot thereof
def TimeEst(time_list, total, output_dir_path, source_dict, kwargs_dict):

    # Add current timing to list of timings, and pass to time estimation function to get predicted completion time
    time_list.append(time.time())
    time_est = ChrisFuncs.TimeEst(time_list, total, plot=True)
    time_remaining = time_est[0]

    # Write estimated completion time to text file
    time_file = open( os.path.join(output_dir_path,'Estimated_Completion_Time.txt'), 'w')
    time_file.write(time_remaining)
    time_file.close()

    # Make plot showing timings so far, and predicted time remaining
    time_fig = time_est[1]
    time_fig.savefig( os.path.join(output_dir_path,'Estimated_Completion_Time.png'), dpi=150  )
    time_fig.clf()
    plt.close('all')

    # If vorbose, report estimated time until completion to user
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] CAAPR estimated completion at: '+time_remaining+'.'



