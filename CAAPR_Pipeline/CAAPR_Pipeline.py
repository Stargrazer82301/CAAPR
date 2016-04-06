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

            # Update current row of bands table to reflect the path of the freshly-made cutout
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
            del(pool)
            pool = mp.Pool(processes=kwargs_dict['n_proc'])
            for band in bands_dict_keys:
                photom_output_list.append( pool.apply_async( CAAPR_Photom.PipelinePhotom, args=(source_dict, bands_dict[band], kwargs_dict) ) )
            pool.close()
            pool.join()
            pdb.set_trace()
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
    [os.remove(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',processed_map)) for processed_map in os.listdir(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps')) if '.fits' in processed_map]
    if os.path.exists(os.path.join(kwargs_dict['temp_dir_path'],'Cutouts',source_dict['name'])):
        shutil.rmtree(os.path.join(kwargs_dict['temp_dir_path'],'Cutouts',source_dict['name']))
    gc.collect()





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
    poly_full = congrid.congrid(poly_fit, (pod['cutout'].shape[0], pod['cutout'].shape[1]), minusone=True)



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




"""
# Define function that runs a map through AstroMagic
from pts.magic import ImageImporter, Extractor
from astropy.units import Unit
from pts.core.tools import logging
def AstroMagic(pod):
    in_fitspath = pod['in_fitspath']
    temp_dir_path = pod['temp_dir_path']
    band_dict = pod['band_dict']
    source_dict = pod['source_dict']



    # Warn user that AstroMagic not yet in place, so placeholder will be used.
    placeholder = True
    if placeholder:
        print 'ASTROMAGIC PROCESSING NOT CURRENTLY ENABLED, PENDING FULL ASTROMAGIC RELEASE. RETURNING UNALTERED MAP AS PLACEHOLDER.'
        am_output = astropy.io.fits.getdata(in_fitspath)

    # If AstroMagic is in place, go ahead (but with output supressed)
    elif not placeholder:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sys.stdout = open(os.devnull, "w")

                # The path to the log file (if you want to have this)
                if os.path.exists(os.path.join(temp_dir_path,band_dict['band_name']))==False:
                    os.mkdir(os.path.join(temp_dir_path,band_dict['band_name']))
                logfile_path = os.path.join(temp_dir_path,band_dict['band_name'],"astromagic_log.txt")

                # Setup the logger
                level = "SUCCESS" # "DEBUG" or "INFO", "WARNING", "ERROR", "SUCCESS"
                logging.setup_log(level=level, path=logfile_path)

                # The path to the image (relative or absolute)
                image_path = in_fitspath

                # If you know the FWHM of the image beforehand, you can specify it for the ImageImporter, and the resulting Image object will
                # have a fwhm attribute. This will then be recognized by the Extractor and it will skip the fitting procedure.
                fwhm = band_dict['beam_arcsec'] * Unit("arcsec") # this should be an Astropy quantity

                # If there are bad pixels in the image, you can create a region file with shapes that cover these pixels.
                # Important is that this region is in pixel coordinates.
                # Examples of what I see as bad pixels is very negative values at the edges, strange artefacts in the image such as
                # stripes ... Pixels with value NaN in the image are not a problem, you don't have to indicate those. These are automatically
                # ignored during the subtraction because the ImageImporter creates a mask from these pixels.
                # Sometimes, saturated stars have a few pixels that are nan. In this case, we certainly don't want to ignore these pixels
                # because we want to remove the star and its diffraction spikes. So, we want to remove these little blobs of nan from the nan_mask,
                # and interpolate the image there ... So here we seperate the nan_mask into a new mask without the little blobs, and a mask consisting
                # of only these blobs. This is all done by the ImageImporter class. If you don't want to set a region of bad pixels manually, just
                # leave the bad_region_path to None.
                bad_region_path = None

                # Import the image
                importer = ImageImporter()
                importer.run(image_path, bad_region_path=bad_region_path, fwhm=fwhm)

                # Create an Extractor instance
                extractor = Extractor()

                # Set configuration options for the extractor

                # Input and output
                #extractor.config.input_path = "test_in" # If there is input (not the image itself), this is where Extractor can find it (you don't need this)
                extractor.config.output_path = os.path.join(temp_dir_path,band_dict['band_name']) # This is where all output goes

                # Here are some options you don't have to enable, but can provide useful output
                # These options are all disabled by default
                extractor.config.write_regions = False
                extractor.config.write_mask = False
                extractor.config.writing.mask_path = os.path.join(temp_dir_path,band_dict['band_name'],"mask.fits")
                extractor.config.write_galactic_catalog = False
                extractor.config.writing.galactic_catalog_path = os.path.join(temp_dir_path,band_dict['band_name'],"galaxies.cat")
                extractor.config.write_stellar_catalog = False
                extractor.config.writing.stellar_catalog_path = os.path.join(temp_dir_path,band_dict['band_name'],"stars.cat")
                extractor.config.write_segments = False

                # This is important for the end of the script:
                extractor.config.writing.result_path = os.path.join( temp_dir_path,  source_dict['name']+'_'+band_dict['band_name']+'_StarSub.fits')

                # You can change the threshold for removing stars from the catalog with the option below
                extractor.config.stars.removal.method = ["interpolation", "interpolation", None]
                # The first element is used for stars which could be fitted to a PSF, the second one is for stars that could not be fitted but for which a peak was detected,
                # and the third is for neither (if you replace None by "interpolation" it will basically interpolate over all positions in the 2MASS point sources catalog).
                # For positions in this catalog where fitting failed, a default FWHM is calculated based on the stars that could be fitted.
                # If the FWHM was specified for the image (see above), then this default FWHM is exactly this specified FWHM. Note that if you specify the FWHM, then not a single
                # star will be fitted to a PSF, so the first option here is ignored. I also want to note that I wanted to allow an option "model" for the first entry, but this
                # is not working properly yet so just set each of these 3 entries to either "interpolation" or None.

                # As a said in one of my previous e-mails, the saturation detection is much too sensitive now.
                #extractor.config.stars.find_saturation = False # this disables the saturation detection completely

                # Here are some parameters I'm playing with to try to improve the saturation detection/removal
                extractor.config.stars.saturation.sigma_level = 5. # the default is 1.2 now ...
                extractor.config.stars.saturation.only_brightest = True
                extractor.config.stars.saturation.brightest_method = "percentage"
                extractor.config.stars.saturation.brightest_level = 10.0

                # Run the extraction (galaxy extraction, star extraction (based on 2MASS point sources) and 'all other sources' extraction (by image segmentation) are all performed here)
                extractor.run(importer.image)

                # Save the result (with the original header if this is desirable)
                extractor.write_result(importer.image.original_header)

                # Grab processed file and pass to pod
                am_output = astropy.io.fits.getdata( os.path.join( temp_dir_path,  source_dict['name']+'_'+band_dict['band_name']+'_StarSub.fits') )

                 # Restore console output, and record map to pod
                sys.stdout = sys.__stdout__
                pod['cutout'] = am_output

        # If AstroMagic failed, report to user, and keep regular image in pod
        except:
            sys.stdout = sys.__stdout__
            print '['+pod['id']+'] Astromagic processing failed; retaining standard image.'

    # Return pod
    return pod
"""