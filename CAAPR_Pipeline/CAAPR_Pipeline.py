# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
sys.path.append( str( os.path.join( os.path.split( os.path.dirname(os.path.abspath(__file__)) )[0], 'CAAPR.CAAPR_AstroMagic', 'PTS') ) )
import gc
import pdb
import time
import re
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
import CAAPR





# The main pipeline; the cutout-production, aperture-fitting, and actual photometry parts of the CAAPR process are called in here, as sub-pipelines
def PipelineMain(source_dict, bands_dict, kwargs_dict):



    # Start timer, and check that the user has actually asked CAAPR to do something; if they haven't asked CAAPR to do anything at all, tell them that they're being a bit odd!
    source_start = time.time()
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Processing target '+source_dict['name']+'.'
    if kwargs_dict['fit_apertures']==False and kwargs_dict['do_photom']==False:
        print '['+source_dict['name']+'] So you don\'t want aperture fitting, nor do you want actual photometry to happen? Erm, okay.'



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

        # In standard operation, process multiple sources in parallel
        aperture_start = time.time()
        aperture_output_list = []
        if kwargs_dict['parallel']==True:
            bands_dict_keys = bands_dict.keys()
            random.shuffle(bands_dict_keys)
            pool = mp.Pool(processes=kwargs_dict['n_proc'])
            for band in bands_dict_keys:
                aperture_output_list.append( pool.apply_async( CAAPR.CAAPR_Aperture.PipelineAperture, args=(source_dict, bands_dict[band], kwargs_dict) ) )
            pool.close()
            pool.join()
            del(pool)
            aperture_list = [output.get() for output in aperture_output_list if output.successful()==True]
            aperture_list = [aperture for aperture in aperture_list if aperture!=None]

        # If parallelisation is disabled, process sources one-at-a-time
        elif kwargs_dict['parallel']==False:
            for band in bands_dict.keys():
                aperture_output_list.append( CAAPR.CAAPR_Aperture.PipelineAperture(source_dict, bands_dict[band], kwargs_dict) )
                aperture_list = [output for output in aperture_output_list if output!=None]

        # Combine all fitted apertures to produce amalgam aperture
        aperture_combined = CAAPR.CAAPR_Aperture.CombineAperture(aperture_list, source_dict, kwargs_dict)

        # Record aperture properties to file
        CAAPR.CAAPR_IO.RecordAperture(aperture_combined, source_dict, kwargs_dict)

        # Create grid of thumbnail images
        CAAPR.CAAPR_IO.ApertureThumbGrid(source_dict, bands_dict, kwargs_dict, aperture_list, aperture_combined)

        # Report time taken to fit apertures, and tidy up garbage
        gc.collect()
        if kwargs_dict['verbose']: print '['+source_dict['name']+'] Time taken performing aperture fitting: '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(time.time()-aperture_start,4))+' seconds.'




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
                photom_output_list.append( pool.apply_async( CAAPR.CAAPR_Photom.PipelinePhotom, args=(source_dict, bands_dict[band], kwargs_dict) ) )
            pool.close()
            pool.join()
            photom_list = [output.get() for output in photom_output_list if output.successful()==True]
            photom_list = [photom for photom in photom_list if photom!=None]

        # If parallelisation is disabled, process sources one-at-a-time
        elif kwargs_dict['parallel']==False:
            for band in bands_dict.keys():
                photom_output_list.append( CAAPR.CAAPR_Photom.PipelinePhotom(source_dict, bands_dict[band], kwargs_dict) )
                photom_list = [output for output in photom_output_list if output!=None]

        # Record photometry results to file
        CAAPR.CAAPR_IO.RecordPhotom(photom_list, source_dict, kwargs_dict)

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





# Define function that does basic initial handling of band parameters
def BandInitiate(band_dict):



    # Make sure band has content
    if band_dict==None:
        return band_dict

    # Parse band cutout request, converting string to boolean if necessary
    if band_dict['make_cutout']=='True':
        band_dict['make_cutout']=True
    if band_dict['make_cutout']=='False':
        band_dict['make_cutout']=False

    # Reset band directory to inviolate value, to purge any holdovers from previous source
    band_dict['band_dir'] = band_dict['band_dir_inviolate']

    # Return band dict
    return band_dict





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
        if 'band_dir_inviolate' in band_dict.keys():
            band_dict['band_dir'] = band_dict['band_dir_inviolate']
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





# Define function to determine what particular band a given band name refers to
def BandParse(band_name_target):



    # Define dictionary containing lists of possible alternate names for each band
    band_names_dict = {'GALEX_FUV':['GALEX_FUV','FUV','FUV-band','GALEX_f','f','f-band'],
                       'GALEX_NUV':['GALEX_NUV','NUV','NUV-band','GALEX_n','n','n-band'],
                       'SDSS_u':['SDSS_u','u','u-band','SDSS_u-band'],
                       'SDSS_g':['SDSS_g','g','g-band','SDSS_g-band'],
                       'SDSS_r':['SDSS_r','r','r-band','SDSS_r-band'],
                       'SDSS_i':['SDSS_i','i','i-band','SDSS_i-band'],
                       'SDSS_z':['SDSS_z','z','z-band','SDSS_z-band'],
                       'CTIO_U':['CTIO_U','CTIO_U-band'],
                       'CTIO_B':['CTIO_B','CTIO_B-band','B','B-band'],
                       'CTIO_V':['CTIO_V','CTIO_V-band','V','V-band'],
                       'CTIO_R':['CTIO_R','CTIO_R-band'],
                       'CTIO_I':['CTIO_I','CTIO_I-band'],
                       'DSS_B':['DSS_B','DSS1_B','DSSI_B','DSS2_B','DSSII_B','DSS_B-band','DSS1_B-band','DSSI_B-band','DSS2_B-band','DSSII_B-band','DSS_G','DSS1_G','DSSI_G','DSS2_G','DSSII_G','DSS_G-band','DSS1_G-band','DSSI_G-band','DSS2_G-band','DSSII_G-band'],
                       'DSS_R':['DSS_R','DSS1_R','DSSI_R','DSS2_R','DSSII_R','DSS_R-band','DSS1_R-band','DSSI_R-band','DSS2_R-band','DSSII_R-band'],
                       'DSS_I':['DSS_I','DSS1_I','DSSI_I','DSS2_I','DSSII_I','DSS_I-band','DSS1_I-band','DSSI_I-band','DSS2_I-band','DSSII_I-band'],
                       '2MASS_J':['2MASS_J','J','J-band','2MASS_J-band'],
                       '2MASS_H':['2MASS_H','H','H-band','2MASS_H-band'],
                       '2MASS_Ks':['2MASS_Ks','Ks','Ks-band','2MASS_Ks-band','2MASS_K','2MASS_K-band'],
                       'UKIRT_J':['UKIRT_J','UKIRT_J-band','UKIDSS_J','UKIDSS_J-band'],
                       'UKIRT_H':['UKIRT_H','UKIRT_H-band','UKIDSS_H','UKIDSS_H-band'],
                       'UKIRT_K':['UKIRT_K','UKIRT_K-band','K','K-band','UKIDSS_K','UKIDSS_K-band'],
                       'Spitzer_3.6':['Spitzer_3.6','Spitzer_3.6um','Spitzer_3.6mu','Spitzer_IRAC_3.6','Spitzer_IRAC_3.6um','Spitzer_IRAC_3.6mu','Spitzer_IRAC1','Spitzer_I1','IRAC_3.6','IRAC_3.6um','IRAC_3.6mu','IRAC1','I1','Spitzer_IRAC1-band','IRAC1-band','I1-band','3.6','3.6um','3.6mu'],
                       'Spitzer_4.5':['Spitzer_4.5','Spitzer_4.5um','Spitzer_4.5mu','Spitzer_IRAC_4.5','Spitzer_IRAC_4.5um','Spitzer_IRAC_4.5mu','Spitzer_IRAC2','Spitzer_I2','IRAC_4.5','IRAC_4.5um','IRAC_4.5mu','IRAC2','I2','Spitzer_IRAC2-band','IRAC2-band','I2-band','4.5','4.5um','4.5mu'],
                       'Spitzer_5.8':['Spitzer_5.8','Spitzer_5.8um','Spitzer_5.8mu','Spitzer_IRAC_5.8','Spitzer_IRAC_5.8um','Spitzer_IRAC_5.8mu','Spitzer_IRAC3','Spitzer_I3','IRAC_5.8','IRAC_5.8um','IRAC_5.8mu','IRAC3','I3','Spitzer_IRAC3-band','IRAC3-band','I3-band','5.8','5.8um','5.8mu'],
                       'Spitzer_8.0':['Spitzer_8.0','Spitzer_8.0um','Spitzer_8.0mu','Spitzer_IRAC_8.0','Spitzer_IRAC_8.0um','Spitzer_IRAC_8.0mu','Spitzer_IRAC4','Spitzer_I4','IRAC_8.0','IRAC_8.0um','IRAC_8.0mu','IRAC4','I4','Spitzer_IRAC4-band','IRAC4-band','I4-band','8.0','8.0um','8.0mu','Spitzer_8','Spitzer_8m','Spitzer_8mu','Spitzer_IRAC_8','Spitzer_IRAC_8um','Spitzer_IRAC_8mu','IRAC_8','IRAC_8um','IRAC_8mu','8','8um','8mu'],
                       'WISE_3.4':['WISE_3.4','WISE_3.4um','WISE_3.4mu','WISE1','WISE1-band','W1','W1-band','WISE_W1','WISE_W1-band'],
                       'WISE_4.6':['WISE_4.6','WISE_4.6um','WISE_4.6mu','WISE2','WISE2-band','W2','W2-band','WISE_W2','WISE_W2-band']}



    # Loop over alternate band name dictionary entries
    band_altnames_matches = []
    for band_name_key in band_names_dict.keys():
        for band_altname in band_names_dict[band_name_key]:

            # Make band names all-lowercase and alphanumeric-only, for ease of comparison
            band_name_target_comp = re.sub(r'\W+', '', band_name_target).replace('_','').lower()
            band_altname_comp = re.sub(r'\W+', '', band_altname).replace('_','').lower()

            # If target and alternate band names match, record
            if band_name_target_comp==band_altname_comp:
                band_altnames_matches.append(band_name_key)

    # If no matches found, or more than one match found, report null output
    if len(band_altnames_matches)==0:
        return None
    elif len(band_altnames_matches)>1:
        return None

    # Else if a good match is found, return it
    elif len(band_altnames_matches)==1:
        return band_altnames_matches[0]




