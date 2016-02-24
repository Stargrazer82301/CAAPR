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
import numpy as np
import scipy.ndimage
import multiprocessing as mp
import astropy.io.fits
sys.path.append( os.path.join( os.path.split( os.path.dirname(os.path.abspath(__file__)) )[0], 'AstroMagic','astromagic') )
sys.path.append( os.path.join( os.path.split( os.path.dirname(os.path.abspath(__file__)) )[0], 'AstroMagic','PTS') )
import congrid
import ChrisFuncs
import astromagic
import pts
from astromagic.core import Image
from astromagic.magic import GalaxyExtractor
from astromagic.magic import StarExtractor
import CAAPR_IO
import CAAPR_Aperture





# The main pipeline; the cutout-production, aperture-fitting, and actual photometry parts of the CAAPR process are called in here, as sub-pipelines
def PipelineMain(source_dict, bands_dict, output_dir_path, temp_dir_path, fit_apertures, aperture_table_path, parallel, n_cores, thumbnails, verbose):



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
            band_cutout_dir = CAAPR_IO.Cutout(source_dict, bands_dict[band], output_dir_path, temp_dir_path)

            # Update current row of bands table to reflect the path of the freshly-made cutout
            bands_dict[band]['band_dir'] = band_cutout_dir



    # If aperture file not provided, commence aperture-fitting sub-pipeline
    if fit_apertures==True:

        # In standard operation, process multiple sources in parallel
        aperture_start = time.time()
        aperture_output_list = []
        if parallel==True:
            bands_dict_keys = bands_dict.keys()
            random.shuffle(bands_dict_keys)
            pool = mp.Pool(processes=n_cores)
            for band in bands_dict_keys:
                aperture_output_list.append( pool.apply_async( CAAPR_Aperture.PipelineAperture, args=(source_dict, bands_dict[band], output_dir_path, temp_dir_path, thumbnails, verbose) ) )
            pool.close()
            pool.join()
            aperture_list = [output.get() for output in aperture_output_list if output.successful()==True]
            aperture_list = [aperture for aperture in aperture_list if aperture!=None]

        # If parallelisation is disabled, process sources one-at-a-time
        elif parallel==False:
            for band in bands_dict.keys():
                aperture_output_list.append( CAAPR_Aperture.PipelineAperture(source_dict, bands_dict[band], output_dir_path, temp_dir_path, thumbnails, verbose) )
                aperture_list = [output for output in aperture_output_list if output!=None]

        # Combine all fitted apertures to produce amalgam aperture
        aperture_combined = CAAPR_Aperture.CombineAperture(aperture_list, source_dict, verbose)
        if verbose: print '['+source_dict['name']+'] Time taken performing aperture fitting: '+str(time.time()-aperture_start)

        # Record aperture properties to file
        aperture_string = str([ source_dict['name'], aperture_combined[0], aperture_combined[1], aperture_combined[2] ])#'name','semimaj_arcsec,axial_ratio,pos_angle\n'
        aperture_string = aperture_string.replace('[','').replace(']','').replace(' ','').replace('\'','')+'\n'
        aperture_table_file = open( aperture_table_path, 'a')
        aperture_table_file.write(aperture_string)
        aperture_table_file.close()

        # Create thumbnail image





# Define function that determines preliminary map values
def MapPrelim(pod, verbose=False):
    band_dict = pod['band_dict']
    source_dict = pod['source_dict']
    verbose = pod['verbose']



    # Check if x & y pixel sizes are meaningfully different. If so, panic; else, treat as same
    pix_size = 3600.0 * pod['in_wcs'].wcs.cdelt
    if float(abs(pix_size.max()))/float(abs(pix_size.min()))>(1+1E-3):
        raise ValueError('The x pixel size if noticably different from the y pixel size.')
    else:
        pod['pix_size'] = float(np.mean(np.abs(pix_size)))

    # Determine source position in cutout in ij coordinates, and size of cutout
    centre_xy = pod['in_wcs'].wcs_world2pix( np.array([[ source_dict['ra'], source_dict['dec'] ]]), 0 )
    pod['centre_i'], pod['centre_j'] = float(centre_xy[0][1]), float(centre_xy[0][0])
    pod['box_rad'] = int( round( float(pod['cutout'].shape[0]) * 0.5 ) )

    # Determine beam size in pixels; if beam size not given, then assume map is Nyquist sampled (ie, 2.355 pixels ber beam)
    if isinstance(band_dict['beam_width'], numbers.Number):
        pod['beam_pix'] = float(band_dict['beam_width']) / pod['pix_size']
    else:
        pod['beam_pix'] = pod['pix_size'] * 2.355

    # Check if current source lies within bounds of map; if not, fai and return)
    if pod['centre_i']<0 or pod['centre_i']>(pod['cutout'].shape)[0] or pod['centre_j']<0 or pod['centre_j']>(pod['cutout'].shape)[1]:
        pod['within_bounds'] = False
        if verbose: print '['+pod['id']+'] Target not within bounds of map.'
    else:
        pod['within_bounds'] = True

    # Return pod
    return pod





# Define function that runs a map through AstroMagic
def AstroMagic(pod):
    in_fitspath = pod['in_fitspath']
    temp_dir_path = pod['temp_dir_path']
    band_dict = pod['band_dict']
    source_dict = pod['source_dict']



    # Warn user that AstroMagic not yet in place, so placeholder will be used.
    placeholder = True
    if placeholder:
        warnings.warn('ASTROMAGIC PROCESSING NOT CURRENTLY ENABLED, PENDING FULL ASTROMAGIC RELEASE. RETURNING UNALTERED MAP AS PLACEHOLDER.')
        am_output = astropy.io.fits.getdata(in_fitspath)
    elif not placeholder:

        # Load input image
        am_image = Image(in_fitspath)

        # Create GalaxyExtractor and StarExtractor objects
        am_galaxyex = GalaxyExtractor()
        am_starex = StarExtractor()

        # Run the galaxy extraction on the primary image frame
        am_galaxyex.run(am_image.frames.primary)

        # Run the star extraction on the primary image frame, passing the galaxy extractor
        am_starex.run(am_image.frames.primary, am_galaxyex)

        # Save the star-subtracted image as a temporary FITS file, then extract just the data to pass back to the main pipeline in the pod
        am_image.save( os.path.join( temp_dir_path,  source_dict['name']+'_'+band_dict['band_name']+'_StarSub.fits') )
        am_output = astropy.io.fits.getdata( os.path.join( temp_dir_path,  source_dict['name']+'_'+band_dict['band_name']+'_StarSub.fits') )

    # Record processed map to pod, and return
    pod['in_image'] = am_output
    return pod





# Define function that fits and subtracts polynomial background filter from map
def PolySub(pod, poly_order=5, cutoff_sigma=2.0):
    if pod['verbose']: print '['+pod['id']+'] Determining if (and how) background is significnatly variable.'



    # Define Keflavich function to downsample an array
    def Downsample(myarr,factor,estimator=np.nanmean):
        ys,xs = myarr.shape
        crarr = myarr[:ys-(ys % int(factor)),:xs-(xs % int(factor))]
        dsarr = estimator( np.concatenate([[crarr[i::factor,j::factor]
            for i in range(factor)]
            for j in range(factor)]), axis=0)
        return dsarr



    # Determine size of central region to mask (twice the radius of the significant-pixel-enclosing ellipse from CAAPR_Aperture.ApertureShape)
    if 'opt_semimaj' in pod.keys():
        pdb.set_trace()
    else:
        mask_semimaj = 2.0 * pod['semimaj_initial_pix']
        mask_axial_ratio = pod['opt_axial_ratio']
        mask_angle = pod['opt_angle']

    # If image has pixels smaller than some limit, downsample image to improve processing time
    pix_size = pod['pix_size']
    pix_size_limit = 2.0
    if pix_size<pix_size_limit:
        downsample_factor = int(np.ceil(pix_size_limit/pix_size))
    else:
        downsample_factor = 1
    image_ds = Downsample(pod['cutout'], downsample_factor)

    # Downsample related values accordingly
    mask_semimaj = mask_semimaj / downsample_factor
    centre_i = int(round(float((0.5*pod['centre_i'])-1.0)))
    centre_j = int(round(float((0.5*pod['centre_j'])-1.0)))



    # Find cutoff for excluding bright pixels by sigma-clipping map
    clip_value = ChrisFuncs.SigmaClip(image_ds, tolerance=0.01, sigma_thresh=3.0, median=True)
    noise_value = clip_value[0]
    field_value = clip_value[1]
    cutoff = field_value + ( cutoff_sigma * noise_value )



    # Mask all image pixels in masking region around source
    image_masked = image_ds.copy()
    ellipse_mask = ChrisFuncs.Photom.EllipseMask(image_ds, mask_semimaj, mask_axial_ratio, mask_angle, centre_i, centre_j)
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
    poly_full = congrid.congrid(poly_fit, (pod['cutout'].shape[0], pod['cutout'].shape[1]), minusone=True)



    # Do a memory check before continuing
    CAAPR_IO.MemCheck(pod)

    # Establish background variation before application of filter
    clip_in = ChrisFuncs.SigmaClip(pod['cutout'], tolerance=0.005, median=True, sigma_thresh=2.0)
    bg_in = pod['cutout'][ np.where( pod['cutout']<clip_in[1] ) ]
    spread_in = np.mean( np.abs( bg_in - clip_in[1] ) )

    # How much reduction in background variation there was due to application of the filter
    image_sub = pod['cutout'] - poly_full
    clip_sub = ChrisFuncs.SigmaClip(image_sub, tolerance=0.005, median=True, sigma_thresh=2.0)
    bg_sub = image_sub[ np.where( image_sub<clip_sub[1] ) ]
    spread_sub = np.mean( np.abs( bg_sub - clip_sub[1] ) )
    spread_diff = spread_in / spread_sub

    # If the filter made significant difference, apply to image and return it; otherwise, just return the unaltered map
    if spread_diff>1.1:
        if pod['verbose']: print '['+pod['id']+'] Background is significnatly variable; removing polynomial background fit, then re-determining source shape.'
        pod['cutout'] = image_sub
        pod['sky_poly'] = poly_model
    else:
        if pod['verbose']: print '['+pod['id']+'] Background is not significnatly variable; leaving image unaltered.'
        pod['sky_poly'] = False
    return pod



