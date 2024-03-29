# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import time
import math
import random
import re
import copy
import warnings
import multiprocessing as mp
import numpy as np
import scipy.optimize
import scipy.ndimage.measurements
import scipy.stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.io.fits
import astropy.wcs
import astropy.convolution
import astropy.modeling
import skimage.restoration
import lmfit
import skimage
import ChrisFuncs
import ChrisFuncs.Photom
import ChrisFuncs.FromGitHub
import CAAPR
plt.ioff()





# The aperture-fitting sub-pipeline
def SubpipelinePhotom(source_dict, band_dict, kwargs_dict):
    source_id = source_dict['name']+'_'+band_dict['band_name']



    # Carry out small random wait, to stop RAM checks from syncing up
    time.sleep(5.0*np.random.rand())



    # Perform initial checks of target file type and location; return if not present
    in_fitspath, file_found = CAAPR.CAAPR_Pipeline.FilePrelim(source_dict, band_dict, kwargs_dict)
    if file_found==False:
        output_dict = {'band_name':band_dict['band_name'], 'ap_sum':np.NaN, 'ap_error':np.NaN}
        return output_dict



    # Create the pod (Photometry Organisation Dictionary), which will read in the FITS file, and bundle all the photometry data for this source & band into one dictionary to be passed between functions
    pod = CAAPR.CAAPR_Pipeline.PodInitiate(in_fitspath, source_dict, band_dict, kwargs_dict)



    # Run pod through preliminary processing, to determine initial quantities; if target not within bounds of map, end processing here
    pod = CAAPR.CAAPR_Pipeline.MapPrelim(pod, source_dict, band_dict)
    if pod['within_bounds']==False:
        output_dict = {'band_name':band_dict['band_name'], 'ap_sum':np.NaN, 'ap_error':np.NaN}
        return output_dict
    CAAPR.CAAPR_IO.MemCheck(pod)



    # Read in aperture file, extract aperture for current source, and apply aperture corrections for beam size
    pod = AperturePrelim(pod, source_dict, band_dict, kwargs_dict)



    # Check if this band is to be excluded from aperture-fitting; if so, return null aperture information
    pod = ExcludePhotom(pod, source_dict, band_dict, kwargs_dict)
    if pod['band_exclude']==True:
        return pod['null_output_dict']



    # If star-removal is required, run pod through AstroMagic
    pod = CAAPR.CAAPR_AstroMagic.Magic(pod, source_dict, band_dict, kwargs_dict)



    # Run pod through function that removes large-scale sky using a 2-dimensional polynomial filter, with source aperture masked
    pod = CAAPR.CAAPR_Pipeline.PolySub(pod, pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'], instant_quit=max([not kwargs_dict['polysub'],pod['band_exclude']]))

    # Processed images as final product: do extinction correction, then save to output_dir
    if kwargs_dict['save_images']:
        # Use IRSA dust extinction service to correction images for extinction
        pod = ExtCorrrct(pod, source_dict, band_dict, kwargs_dict)

        astropy.io.fits.writeto(os.path.join(kwargs_dict['output_dir_path'],'Processed_Maps',source_id+'.fits'), pod['cutout'], header=pod['in_header'], overwrite=True)

    # If no photometry is done, we can return early
    if not kwargs_dict['do_photom']:
        return

    # Run pod through function that (finally) performs the actual photometry
    pod = Photom(pod, band_dict)



    # If photometry recorded null flux, skip determination of aperture noise, recording null value for this too, instead
    if np.isnan(pod['ap_sum'])==True:
        pod['ap_error'] = np.NaN
    else:
        CAAPR.CAAPR_IO.MemCheck(pod)



        # Attempt to determine aperture noise the preferred way, using full-size randomly-placed apertures
        if kwargs_dict['verbose']: print '['+source_id+'] Estimating aperture noise using full-size randomly-placed sky apertures.'
        ap_noise_dict = ApNoise(pod['cutout'].copy(), source_dict, band_dict, kwargs_dict, pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'], downsample=int(band_dict['downsample_factor']))
        if ap_noise_dict['fail']==False:
            if kwargs_dict['verbose']: print '['+source_id+'] Aperture noise successfully estimated using full-size randomly-placed sky apertures.'
            pod['ap_noise'] = ap_noise_dict['ap_noise']



        # If full-size apertures were unable to produce a valid aperture noise estimate, commence aperture extrapolation
        else:
            sky_success_counter = ap_noise_dict['sky_success_counter']
            if kwargs_dict['verbose']: print '['+source_id+'] Unable to estiamte aperture noise using full-size randomly-placed sky apertures (only '+str(int(sky_success_counter))+' could be placed); switching to aperture extrapolation.'
            CAAPR.CAAPR_IO.MemCheck(pod)
            ap_noise_dict = ApNoiseExtrap(pod['cutout'].copy(), source_dict, band_dict, kwargs_dict, pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'])
            pod['ap_noise'] = ap_noise_dict['ap_noise']



        # If even aperture extrapolation is unable to produce an aperture noise estimate, report null value
        if ap_noise_dict['fail']==True:
            pod['ap_noise'] = np.NaN
        else:
            if pod['verbose']: print '['+pod['id']+'] Final aperture noise is '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(pod['ap_noise'],4))+' (in map units).'



        # Calculate calibration uncertainty
        calib_err = abs( pod['ap_sum'] * band_dict['calib_error'] )
        if pod['verbose']:print '['+pod['id']+'] Calibration uncertainty is '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(calib_err,4))+' (in map units).'

        # Calculate final uncertainty
        pod['ap_error'] = ( pod['ap_noise']**2.0 + calib_err**2.0 )**0.5
        if np.isnan(pod['ap_noise'])==True:
            pod['ap_error'] = -1.0 * pod['ap_sum'] * band_dict['calib_error']



        # Run pod through function that performs aperture correction
        pod = ApCorrect(pod, source_dict, band_dict, kwargs_dict)



        # Use IRSA dust extinction service to correction fluxes for extinction
        pod = ExtCorrrct(pod, source_dict, band_dict, kwargs_dict)



    # If thumbnail images have been requested, save a copy of the current image (ie, with any star and/or background subtaction)
    if kwargs_dict['thumbnails']==True:
        astropy.io.fits.writeto(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_id+'.fits'), pod['cutout'], header=pod['in_header'], overwrite=True)



    # Now return final photometry informaton to main pipeline, and clean up garbage
    if np.isnan(pod['ap_sum'])==False:
        if pod['verbose']:print '['+pod['id']+'] Final source flux and uncertainty is '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(pod['ap_sum'],4))+' +/- '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(pod['ap_error'],4))+' (in map units).'
    output_dict = {'band_name':band_dict['band_name'], 'ap_sum':pod['ap_sum'], 'ap_error':pod['ap_error']}
    gc.collect()
    del(pod)
    return output_dict





# Define functiont that reads in aperture file, identifies aperture for current source, and applies aperture corrections for beam size
def AperturePrelim(pod, source_dict, band_dict, kwargs_dict):



    # Read in aperture file
    aperture_table = np.genfromtxt(kwargs_dict['aperture_table_path'], delimiter=',', names=True, dtype=None)
    aperture_index = np.where( aperture_table['name']==source_dict['name'] )
    if aperture_index[0].shape[0]>1:
        raise Exception('Aperture value caontains more than one entry for current galaxy')
    else:
        aperture_index = aperture_index[0][0]

    # Extract aperture corresponding to current source, dealing with special case where aperture file contains only one source
    if aperture_table['semimaj_arcsec'].shape==() and np.where( aperture_table['name']==source_dict['name'] )[0][0]==0:
        opt_semimaj_arcsec = aperture_table['semimaj_arcsec']
        opt_axial_ratio = aperture_table['axial_ratio']
        opt_angle = aperture_table['pos_angle']
    else:
        opt_semimaj_arcsec = aperture_table['semimaj_arcsec'][aperture_index]
        opt_axial_ratio = aperture_table['axial_ratio'][aperture_index]
        opt_angle = aperture_table['pos_angle'][aperture_index]
    opt_semimin_arcsec = opt_semimaj_arcsec / opt_axial_ratio

    # Adjust aperture to account for beam size of current band, and record to pod
    adj_semimaj_arcsec = 0.5 * ( (2.0*opt_semimaj_arcsec)**2.0 + band_dict['beam_arcsec']**2.0 )**0.5
    adj_semimin_arcsec = 0.5 * ( (2.0*opt_semimin_arcsec)**2.0 + band_dict['beam_arcsec']**2.0 )**0.5
    adj_axial_ratio = adj_semimaj_arcsec / adj_semimin_arcsec
    adj_angle = opt_angle

    # Save adjusted values to pod
    pod['adj_semimaj_arcsec'] = adj_semimaj_arcsec
    pod['adj_semimin_arcsec'] = adj_semimin_arcsec
    pod['adj_semimaj_pix'] = adj_semimaj_arcsec / pod['pix_arcsec']
    pod['adj_semimin_pix'] = adj_semimin_arcsec / pod['pix_arcsec']
    pod['adj_axial_ratio'] = adj_axial_ratio
    pod['adj_angle'] = adj_angle

    # Return pod
    return pod





# Define actual function that actually performs actual photometry
def Photom(pod, band_dict):
    if pod['band_exclude']==True:
        return pod
    if pod['verbose']: print '['+pod['id']+'] Performing aperture photometry.'



    # Evaluate pixels in source aperture; with consideration of sub-pixels if requested
    if float(band_dict['subpixel_factor'])==1.0:
        ap_calc = ChrisFuncs.Photom.EllipseSum(pod['cutout'], pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'])
    elif float(band_dict['subpixel_factor'])>1.0:
        ap_calc = ChrisFuncs.Photom.EllipseSumUpscale(pod['cutout'], pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'], upscale=band_dict['subpixel_factor'])

    # Evaluate pixels in background annulus; with consideration of sub-pixels if requested
    bg_inner_semimaj_pix = pod['adj_semimaj_pix'] * band_dict['annulus_inner']
    bg_width = (pod['adj_semimaj_pix'] * band_dict['annulus_outer']) - bg_inner_semimaj_pix
    if float(band_dict['subpixel_factor'])==1.0:
        bg_calc = ChrisFuncs.Photom.AnnulusSum(pod['cutout'], bg_inner_semimaj_pix, bg_width, pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'])
    elif float(band_dict['subpixel_factor'])>1.0:
        bg_calc = ChrisFuncs.Photom.AnnulusSumUpscale(pod['cutout'], bg_inner_semimaj_pix, bg_width, pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'], upscale=band_dict['subpixel_factor'])

    # Check what fraction of aperture an annulus pixels are nan
    nan_fail = False
    if np.where(np.isnan(ap_calc[2]))[0].shape[0]==0:
        ap_nan_frac = 0.0
    else:
        ap_nan_frac = float(np.where(np.isnan(ap_calc[2]))[0].shape[0]) / float(ap_calc[1])
    if np.where(np.isnan(bg_calc[2]))[0].shape[0]==0:
        bg_nan_frac = 0.0
    else:
        bg_nan_frac = float(np.where(np.isnan(bg_calc[2]))[0].shape[0]) / float(bg_calc[1])

    # If more than a givn fraction of the pixels inside the source aperture, or the background annulus, are NaN, record NaN flux, otherwise continue as per normal
    ap_nan_thresh = 0.10
    bg_nan_thresh = 0.80
    if ap_nan_frac>ap_nan_thresh:
        if pod['verbose']: print '['+pod['id']+'] More than '+str(int(100.0*ap_nan_thresh))+'% of pixels in source aperture are NaN; recording null flux.'
        nan_fail = True
    elif bg_nan_frac>bg_nan_thresh:
        if pod['verbose']:print '['+pod['id']+'] More than '+str(int(100.0*bg_nan_thresh))+'% of pixels in background are NaN; recording null flux.'
        nan_fail = True
    if nan_fail==True:
        pod['ap_sum'] = np.NaN
        pod['bg_avg'] = np.NaN
        return pod
    else:

        # Calculate background level, and subtract from flux in source aperture to determine source flux
        bg_clip = ChrisFuncs.SigmaClip(bg_calc[2], median=False, sigma_thresh=3.0)
        bg_avg = bg_clip[1] * float(band_dict['subpixel_factor'])**2.0
        ap_sum = ap_calc[0] - (ap_calc[1] * bg_avg)

        # Save values to pod, and return
        pod['ap_sum'] = ap_sum
        pod['bg_avg'] = bg_avg
        if np.isnan(ap_sum):
            if pod['verbose']:print '['+pod['id']+'] Source flux is NaN.'
        else:
            if pod['verbose']:print '['+pod['id']+'] Source flux is '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(ap_sum,4))+' (in map units).'
        return pod





# Define funcion that attempts to estimate aperture noise using randomly-positioned sky apertures of given dimensions
def ApNoise(cutout, source_dict, band_dict, kwargs_dict, adj_semimaj_pix, adj_axial_ratio, adj_angle, centre_i, centre_j, mini=False, downsample=False):
    source_id = source_dict['name']+'_'+band_dict['band_name']
    ap_debug = False



    # Handle downsampling input
    if downsample!=False:
        ds_request = int(downsample)
    else:
        ds_request = 1

    # Standard downsampling target is for aperture diameter to correspod to 200 pixels
    ds_target = int( np.round( float(adj_semimaj_pix) ) / 100.0 )
    ds_factor = max([ ds_request, ds_target ])

    # If mini-apertures are being used, ensure downsampling isn't agressive (ie, apertures would be sub-Nyquist sampled)
    if mini!=False:
        ds_mini = int( np.round( float(mini-1.0) ) / 2.355 )
        if ds_mini>=2:
            ds_factor = min([ ds_factor, ds_mini ])
        else:
            ds_factor = 1

    # If downsampling is makes sense, apply
    if ds_factor>=2:
        if ap_debug: print 'Setup: Downsampling map by factor of '+str(ds_factor)
        cutout = skimage.measure.block_reduce(cutout, block_size=(int(ds_factor),int(ds_factor)), func=np.mean, cval=np.NaN)
        cutout *= float(ds_factor)*float(ds_factor)
        centre_i /= float(ds_factor)
        centre_j /= float(ds_factor)
        adj_semimaj_pix /= float(ds_factor)
        if mini!=False:
            mini /= float(ds_factor)



    # Handle input variables if mini-apertures are required
    if mini!=False:
        if ap_debug: print 'Setup: Preparing inputs for mini-apertures'
        if isinstance(mini, float) or isinstance(mini, int):
            mini = float(mini)
            adj_semimaj_pix_full = adj_semimaj_pix
            adj_semimaj_pix = mini
        else:
            pdb.set_trace()
    else:
        adj_semimaj_pix_full = adj_semimaj_pix
    adj_semimin_pix_full = adj_semimaj_pix_full / adj_axial_ratio



    # Define charactaristics of circular aperture with same area as elliptical source aperture
    ap_area = np.pi * adj_semimaj_pix * ( adj_semimaj_pix / adj_axial_ratio )
    sky_ap_rad_pix = ( ap_area / np.pi )**0.5
    sky_border = int( sky_ap_rad_pix + 1.0 )#int( ( band_dict['annulus_outer'] * sky_ap_rad_pix ) + 1 )
    adj_semimin_pix = adj_semimaj_pix / adj_axial_ratio

    # Creating mask maps to describe no-go regions
    if ap_debug: print 'Setup: Creating mask maps'
    prior_mask = np.zeros(cutout.shape)
    exclude_mask = ChrisFuncs.Photom.EllipseMask(cutout, adj_semimaj_pix_full, adj_axial_ratio, adj_angle, centre_i, centre_j)
    flag_mask = np.zeros(cutout.shape)
    attempt_mask = np.zeros(cutout.shape)

    # Set pixels in source aperture to all have NaN pixels, so they don't get sampled by sky annuli
    cutout_inviolate = cutout.copy()
    cutout[ np.where(ChrisFuncs.Photom.EllipseMask(cutout, adj_semimaj_pix_full, adj_axial_ratio, adj_angle, centre_i, centre_j)==1) ] = np.NaN



    # Define how many random aperture are desired/required/permitted
    sky_success_target = 100
    sky_success_min = 50
    sky_gen_max = 250

    # Generate random polar coordinates to draw from
    if ap_debug: print 'Setup: Generating pool of random polar coordinates'
    random_size = sky_success_target * sky_gen_max * 10
    random_failed = []
    random_theta_list = 360.0 * np.random.rand(random_size)
    random_r_list = adj_semimin_pix + np.abs(np.random.normal(loc=0.0, scale=5.0*adj_semimaj_pix_full, size=random_size))

    # Locate contiguous map regions
    if ap_debug: print 'Pruning: Locating contiguous coverage regions'
    cont_binary = np.zeros(cutout.shape)
    cont_binary[ np.where( np.isnan(cutout_inviolate)==False ) ] = 1
    cont_structure = np.array([[1,1,1], [1,1,1], [1,1,1]])
    cont_label = scipy.ndimage.measurements.label(cont_binary, structure=cont_structure)[0]
    cont_search_mask = ChrisFuncs.EllipseMask(cont_label, 3.0, 1.0, 0.0, centre_i, centre_j)
    cont_search_values = cont_label[ np.where( cont_search_mask==1 ) ]

    # Identify contiguous region associated with target source
    if ap_debug: print 'Pruning: Identifying coverage region associated with source'
    if np.where(cont_search_values>0)[0].shape[0]==0:
        cont_target = 0
    else:
        cont_target = scipy.stats.mode(cont_search_values[np.where(cont_search_values>0)])[0][0]
    cont_where_bad = np.where( cont_label!=cont_target )

    # Remove random coordinates that are more distant than most distant part of the coverage region the target source lies in
    if ap_debug: print 'Pruning: Removing random coords that definately lie outside coverage region'
    cont_size_i, cont_size_j = cutout.shape
    cont_range_i = np.arange(cont_size_i) - centre_i
    cont_range_j = np.arange(cont_size_j) - centre_j
    cont_coord_i, cont_coord_j = np.meshgrid(cont_range_j, cont_range_i) # Yes, i and j are supposed to be this way around inside meshgrid (for some reason)
    cont_coord_i[cont_where_bad] = np.NaN
    cont_coord_j[cont_where_bad] = np.NaN
    cont_dist = np.sqrt(cont_coord_i**2 + cont_coord_j**2)
    cont_dist_max = np.nanmax(cont_dist)
    random_r_coverage = np.where( random_r_list < (cont_dist_max-sky_border) )
    random_r_list = random_r_list[random_r_coverage]
    random_theta_list = random_theta_list[random_r_coverage]

    # Convert random polar coordinates into cartesian coordinates
    if ap_debug: print 'Pruning: Converting random polar coords to cartesian coords, and removing those that lie beyond map border'
    random_i_list = centre_i + ( random_r_list * np.cos(np.radians(random_theta_list)) )#np.random.normal(loc=centre_i, scale=2.0*sky_ap_rad_pix)
    random_j_list = centre_j + ( random_r_list * np.sin(np.radians(random_theta_list)) )#np.random.normal(loc=centre_j, scale=2.0*sky_ap_rad_pix)

    # Remove random coodinates that fall fall beyond border in i-coords
    random_not_i_border = np.where( (random_i_list>sky_border) & (random_i_list<(cutout.shape[0]-sky_border)) )
    random_i_list = random_i_list[random_not_i_border]
    random_j_list = random_j_list[random_not_i_border]

    # Remove random coodinates that fall fall beyond border in j-coords
    random_not_j_border = np.where( (random_j_list>sky_border) & (random_j_list<(cutout.shape[1]-sky_border)) )
    random_i_list = random_i_list[random_not_j_border]
    random_j_list = random_j_list[random_not_j_border]

    # Remove random coordinates that intersect source
    if ap_debug: print 'Pruning: Removing random coords that intersect source'
    random_not_source = np.where( np.sqrt( (np.abs(centre_i-random_i_list))**2.0 + (np.abs(centre_j-random_j_list))**2.0 ) > adj_semimin_pix_full )
    #random_not_source = np.where( (abs(centre_i-random_i_list)>adj_semimin_pix_full) & (abs(centre_j-random_j_list)>adj_semimin_pix_full) )
    random_i_list = random_i_list[random_not_source]
    random_j_list = random_j_list[random_not_source]

    # Remove random coordinates that correspond to NaN pixels
    if ap_debug: print 'Pruning: Removing random coords that correspond to NaN pixels'
    random_i_list_pix = np.round(random_i_list).astype(int)
    random_j_list_pix = np.round(random_j_list).astype(int)
    random_ij_values = cutout[(random_i_list_pix,random_j_list_pix)]
    random_ij_pix_good = np.where(np.isnan(random_ij_values)==False)
    random_i_list = random_i_list[random_ij_pix_good]
    random_j_list = random_j_list[random_ij_pix_good]

    # If none of the apertures are suitable, immediately report failure
    if random_i_list.shape[0]==0:
        if ap_debug: print 'Status: Pruning removed all generated coordinates'
        ap_noise_dict = {'fail':True, 'prior_mask':prior_mask, 'flag_mask':flag_mask, 'sky_success_counter':0}
        cutout = cutout_inviolate
        return ap_noise_dict

    # Commence creation of random sky apertures
    sky_success_counter = 0
    sky_sum_list = []
    sky_total_fail = False
    while True:

        # Repeatedly generate random sky apertures, until an acceptable aperture is generated
        sky_gen_counter = 0
        sky_gen_fail = False
        while True:
            sky_gen_counter += 1

            # If more than a given number of unsuccessful sky apertures have been generated in a row, call it a day
            if sky_gen_counter>sky_gen_max:
                sky_gen_fail = True
                if ap_debug: print 'Status: Unable to generate suitable random sky aperture after '+str(sky_gen_max)+' attempts'
                break



            # Select random coordinate set for this iteration; if no un-used random coordinates can be found, reject
            #if ap_debug: print 'Checking: Selecting random coordinates'
            random_accept = False
            random_reject_count = 0
            while random_accept==False:
                random_index = int( np.floor( np.random.rand() * float(random_i_list.shape[0]) ) )
                if random_index not in random_failed:
                    random_accept = True
                else:
                    random_reject_count += 1
                if random_reject_count>(10*random_i_list.shape[0]):
                    break
            if random_reject_count>(10*random_i_list.shape[0]):
                if ap_debug: print 'Rejection: Unable to find un-used random coodinates'
                continue
            random_i = random_i_list[random_index]
            random_j = random_j_list[random_index]
            if ap_debug:
                ap_mask = ChrisFuncs.Photom.EllipseMask(cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
                attempt_mask[ np.where(ap_mask==1) ] = sky_success_counter
                print 'Aperture: '+str(sky_success_counter+1)+';   Generation: '+str(sky_gen_counter)+';   Pix Coords: ['+str(random_i)+','+str(random_j)+']'

            # Do sophisticated check that generated sky aperture does not intersect source; if it does, reject
            #if ap_debug: print 'Checking: Determining whether aperture intersects source (sophisticated check)'
            exclude_sum = ChrisFuncs.Photom.EllipseSum(exclude_mask, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)[0]
            if exclude_sum>0:
                if ap_debug: print 'Rejection: Aperture intersects source (according to sophisticated check)'
                random_failed.append(random_index)
                continue

            # Do basic chrck that the majority of the pixels in the generated sky aperture have not already been sampled by previous sky apertures; they have, reject
            #if ap_debug: print 'Checking: Determining if aperture over-sampled (basic check)'
            prior_calc = ChrisFuncs.Photom.EllipseSum(prior_mask, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
            prior_calc[2][ np.where(prior_calc[2]>=1.0) ] = 1.0
            prior_frac = np.sum(prior_calc[2]) / float(prior_calc[1])
            if prior_frac>0.5:
                if ap_debug: print 'Rejection: Aperture over-sampled (according to basic check)'
                random_failed.append(random_index)
                continue

            # Do sophisticated check that the majority of the pixels in the generated sky aperture have not already been sampled by previous sky apertures; they have, reject
            #if ap_debug: print 'Checking: Determinging if aperture oversampled (sophisticated check)'
            ap_mask_check = ChrisFuncs.Photom.EllipseMask(cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
            flag_mask_check = flag_mask.copy()
            flag_mask_check[np.where(ap_mask_check==1)] = int(2.0**(sky_success_counter+1.0))
            flag_tallies = np.array([ np.where(flag_mask_check==flag)[0].shape[0] for flag in (2.0**np.arange(0.0,sky_success_counter+2.0)).tolist() ])
            flag_check = np.where(flag_tallies<(0.5*ap_area))[0].shape[0]
            if flag_check>1:
                if ap_debug: print 'Rejection: Aperture over-sampled (according to sophisticated check)'
                random_failed.append(random_index)
                continue

            # Evaluate pixels in sky aperture
            #if ap_debug: print 'Checking: Evaluating pixels in sky aperture'
            ap_calc = ChrisFuncs.Photom.EllipseSum(cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)

            # Evaluate pixels in sky annulus
            #if ap_debug: print 'Checking: Evaluating pixels in sky anmnulus'
            bg_inner_semimaj_pix = adj_semimaj_pix * band_dict['annulus_inner']
            bg_width = (adj_semimaj_pix * band_dict['annulus_outer']) - bg_inner_semimaj_pix
            bg_calc = ChrisFuncs.Photom.AnnulusSum(cutout, bg_inner_semimaj_pix, bg_width, 1.0, 0.0, random_i, random_j)

            # Check if more than a given fraction of the pixels inside the source aperture are NaN; if so, reject
            if ap_calc[3][0].shape[0]==0 or ap_calc[1]==0:
                ap_nan_frac = 0.0
            if ap_calc[1]==0:
                ap_nan_frac = 1.0
            else:
                ap_nan_frac = float(ap_calc[3][0].shape[0]) / float(ap_calc[1]+float(ap_calc[3][0].shape[0]))
            ap_nan_thresh = 0.10
            if ap_nan_frac>ap_nan_thresh:
                if ap_debug: print 'Rejection: Aperture contains too many NaNs'
                random_failed.append(random_index)
                continue

            # Check if more than a given fraction of the pixels inside the sky annulus are NaN; if so, reject
            if bg_calc[3][0].shape[0]==0:
                bg_nan_frac = 0.0
            if bg_calc[1]==0:
                bg_nan_frac = 1.0
            else:
                bg_nan_frac = float(bg_calc[3][0].shape[0]) / float(bg_calc[1]+bg_calc[3][0].shape[0])
            bg_nan_thresh = 0.80
            if bg_nan_frac>bg_nan_thresh:
                if ap_debug: print 'Rejection: Annulus contains too many NaNs'
                random_failed.append(random_index)
                continue



            # If coords have not been rejected for any reason, accept them and proceed
            else:
                sky_success_counter += 1
                break



        # If no suitable sky aperture could be generated on this iteration, decide how to proceed, based on how many had been successfully generated already
        if sky_gen_fail:
            if sky_success_counter<sky_success_min:
                sky_total_fail = True
                break
            else:
                if ap_debug: print 'Status: However, sufficient number of successful random apertures ('+str(int(sky_success_counter))+') already generated; proceeding'
                break

        # Calculate actual flux in sky aperture, and record
        #if ap_debug: print 'Checking: Performing photometry with random sky aperture and annulus'
        bg_clip = ChrisFuncs.SigmaClip(bg_calc[2], median=False, sigma_thresh=3.0)
        bg_avg = bg_clip[1]
        ap_sum = ap_calc[0] - (ap_calc[1] * bg_avg)
        sky_sum_list.append(ap_sum)
        if np.isnan(ap_sum):
            pdb.set_trace()

        # Add this aperture to the prior mask and flag mask
        if not ap_debug: ap_mask = ChrisFuncs.Photom.EllipseMask(cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
        prior_mask += ap_mask
        flag_mask[np.where(ap_mask==1)] += 2.0**sky_success_counter

        # If target number of sky apertures have been processed, break out of loop
        if sky_success_counter>=sky_success_target:
            if ap_debug: print 'Status: Target number of successful random apertures ('+str(int(sky_success_target))+') generated; proceeding'
            break

    # If total failure was encountered, end process and report now
    if sky_total_fail:
        ap_noise_dict = {'fail':True, 'prior_mask':prior_mask, 'flag_mask':flag_mask, 'sky_success_counter':sky_success_counter}
        cutout = cutout_inviolate
        return ap_noise_dict

    # Otherwise, calculate aperture noise using returned aperture values, and return
    else:
        sky_sum_list = np.array(sky_sum_list)
        ap_noise = ChrisFuncs.SigmaClip(sky_sum_list, tolerance=0.001, median=True, sigma_thresh=3.0)[0]
        ap_noise_dict = {'fail':False, 'ap_noise':ap_noise, 'ap_num':sky_success_counter, 'prior_mask':prior_mask, 'flag_mask':flag_mask}
        #ChrisFuncs.Cutout(prior_mask, '/home/saruman/spx7cjc/DustPedia/Prior.fits')
        if kwargs_dict['verbose']: print '['+source_id+'] Aperture noise from current random apertures is '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(ap_noise,4))+' (in map units).'
        cutout = cutout_inviolate
        return ap_noise_dict





# Define funcion that attempts to estimate aperture noise using randomly-positioned sky apertures of given dimensions
def ApNoiseExtrap(cutout, source_dict, band_dict, kwargs_dict, adj_semimaj_pix, adj_axial_ratio, adj_angle, centre_i, centre_j):
    source_id = source_dict['name']+'_'+band_dict['band_name']



    # Define charactaristics of circular aperture with same area as elliptical source aperture
    ap_area = np.pi * adj_semimaj_pix * ( adj_semimaj_pix / adj_axial_ratio )
    sky_ap_rad_pix = ( ap_area / np.pi )**0.5

    # Generate list of mini-aperture sizes to use, and declare result lists
    mini_ap_rad_base = 1.2
    mini_ap_rad_pix_input = mini_ap_rad_base**np.arange( 1.0, np.ceil( math.log( sky_ap_rad_pix, mini_ap_rad_base ) ) )[::-1]
    min_ap_rad_pix_output = []
    mini_ap_noise_output = []
    mini_ap_num_output = []

    # Loop over radii for mini-apertures
    for mini_ap_rad_pix in mini_ap_rad_pix_input:
        if kwargs_dict['verbose']:print '['+source_id+'] Finding aperture noise for mini-apertures of radius '+str(mini_ap_rad_pix)[:6]+' pixels.'
        mini_ap_noise_dict = ApNoise(cutout.copy(), source_dict, band_dict, kwargs_dict, adj_semimaj_pix, adj_axial_ratio, adj_angle, centre_i, centre_j, mini=mini_ap_rad_pix, downsample=int(band_dict['downsample_factor']))

        # If mini-aperture succeeded, record and proceed; else, call it a day
        if mini_ap_noise_dict['fail']==False:
            min_ap_rad_pix_output.append(mini_ap_rad_pix)
            mini_ap_noise_output.append(mini_ap_noise_dict['ap_noise'])
            mini_ap_num_output.append(mini_ap_noise_dict['ap_num'])
        elif mini_ap_noise_dict['fail']==True:
            if kwargs_dict['verbose']:print '['+source_id+'] Unable to place sufficient number of mini-apertures at this radius.'
        if len(mini_ap_noise_output)>=10:
            break

    # Convert output lists into arrays
    min_ap_rad_pix_output = np.array(min_ap_rad_pix_output)
    mini_ap_noise_output = np.array(mini_ap_noise_output)
    mini_ap_num_output = np.array(mini_ap_num_output).astype(float)



    # If insufficient points to make extrapolation, report failure; else proceed
    if min_ap_rad_pix_output.shape[0]<2:
        ap_noise_dict = {'fail':True, 'ap_noise':np.NaN}
        gc.collect()
        return ap_noise_dict
    else:



        # Calculate log of mini aperture area and noise
        log_mini_ap_area = np.log10(np.pi*min_ap_rad_pix_output**2.0)
        log_mini_ap_noise = np.log10(mini_ap_noise_output)

        # Calculate poisson uncertaity on calculated noise
        mini_ap_noise_err_rel = mini_ap_num_output**0.5 / mini_ap_num_output
        mini_ap_noise_err = np.abs( mini_ap_noise_output * mini_ap_noise_err_rel )

        # Weight points according to distance in log space from true aperture area
        mini_ap_noise_err = mini_ap_noise_err * (1.0 + ( ap_area / 10.0**log_mini_ap_area ) )

        # Translate uncertainties into log space
        log_mini_ap_noise_err = ChrisFuncs.LogError(mini_ap_noise_output, mini_ap_noise_err)

        # Define straight-line function, and fit it to points
        def Line(x,m,c):
            return (m*x)+c
        try:
            line_fit = scipy.optimize.curve_fit(Line, log_mini_ap_area, log_mini_ap_noise, sigma=log_mini_ap_noise_err)
        except:
            ap_noise_dict = {'fail':True, 'ap_noise':np.NaN}
            gc.collect()
            return ap_noise_dict
        line_m, line_c = line_fit[0][0], line_fit[0][1]

        # Determine projected aperture noise value
        log_ap_area = np.log10(ap_area)
        log_ap_noise_proj = Line(log_ap_area, line_m, line_c) #10.0**( ( line_m * np.log10(ap_area) ) + line_c )
        ap_noise_proj = 10.0**(log_ap_noise_proj)

        # Generate points for best-fit line
        ax_y_min = np.floor( np.min([ np.min( log_mini_ap_noise - log_mini_ap_noise_err ), log_ap_noise_proj ]) )
        ax_y_max = np.ceil( np.max([ np.max( log_mini_ap_noise + log_mini_ap_noise_err ), log_ap_noise_proj ]) )
        ax_x_min = np.floor( np.min([ np.min(log_mini_ap_area), log_ap_area ]) )
        ax_x_max = np.ceil( np.max([ np.max(log_ap_area), log_ap_area ]) )
        line_x = np.linspace( ax_x_min, ax_x_max, num=10000 )
        line_y = Line(line_x, line_m, line_c)

        # Set up figure & axes
        fig = plt.figure(figsize=(8,6))
        ax_dims = [0.125, 0.125, 0.825, 0.825]
        ax = fig.add_axes(ax_dims)

        # Plot points and best-fit line
        ax.errorbar(log_mini_ap_area, log_mini_ap_noise, yerr=log_mini_ap_noise_err, ecolor='#4D78C9', elinewidth=1.15, capthick=1.15, marker='x', color='#0080FF', markersize=0.0, markeredgewidth=1.15, linewidth=0)
        ax.scatter(log_mini_ap_area, log_mini_ap_noise, c='#4D78C9', marker='o', s=75, linewidths=0, label='Mini-aperture noise values')
        ax.scatter(np.log10(ap_area), log_ap_noise_proj, c='#C03030', marker='H', s=150, linewidths=0, label='Extrapolated aperture noise')
        ax.plot(line_x, line_y, ls='--', lw=1.0, c='#4D78C9', label='Line of best fit')

        # Format axis limts and labels
        ax.set_xlabel(r'Aperture Area (log$_{10}$ pix)', fontsize=15)
        ax.set_ylabel(r'Aperture Noise (log$_{10}$ map units)', fontsize=15)
        ax.set_xlim(ax_x_min,ax_x_max)
        ax.set_ylim(ax_y_min,ax_y_max)
        for xlabel in ax.get_xticklabels():
            xlabel.set_fontproperties(matplotlib.font_manager.FontProperties(size=15))
        for ylabel in ax.get_yticklabels():
            ylabel.set_fontproperties(matplotlib.font_manager.FontProperties(size=15))

        # Plot axis 1 legend
        ax_handles, ax_labels = ax.get_legend_handles_labels()
        ax1_lgd = ax.legend(ax_handles, ax_labels, loc='best', scatterpoints=1, labelspacing=0.25, borderpad=0)
        ax1_lgd.draw_frame(False)
        plt.setp(plt.gca().get_legend().get_texts(), fontsize='12.5')



        # Save figure, clean up, and report results
        fig.savefig( os.path.join( kwargs_dict['temp_dir_path'], source_id+'_Aperture_Noise_Projection.png' ), dpi=100 )
        gc.collect()
        fig.clear()
        plt.close('all')
        ap_noise_dict = {'fail':False, 'ap_noise':ap_noise_proj}
        return ap_noise_dict





# Define function that uses provided beam profile to aperture-correct photometry
def ApCorrect(pod, source_dict, band_dict, kwargs_dict):



    # If aperture correction not required, immediately return unaltered pod
    if str(band_dict['beam_correction'])=='False':
        return pod
    elif kwargs_dict['verbose']: print '['+pod['id']+'] Determining aperture correction factor due to PSF.'



    # If SNR ratio is <3, do not perform aperture correction
    snr = pod['ap_sum'] / abs(pod['ap_error'])
    if snr<3.0:
        if kwargs_dict['verbose']: print '['+pod['id']+'] Source SNR<3; not applying aperture correction, as reuslt would be unreliable.'
        return pod



    # If no PSF given, assume Airy disc; else extract PSF from provided file
    if str(band_dict['beam_correction'])=='True':
        psf = astropy.convolution.kernels.AiryDisk2DKernel(pod['beam_pix']).array
    else:

        # Read in PSF, and establish pixel size
        psf_in, psf_header = astropy.io.fits.getdata(band_dict['beam_correction'], header=True)
        psf_wcs = astropy.wcs.WCS(psf_header)
        if psf_wcs.wcs.has_cd():
            psf_cdelt = psf_wcs.wcs.cd.max()
        else:
            psf_cdelt = psf_wcs.wcs.cdelt.max()
        psf_cdelt_arcsec = abs( psf_cdelt * 3600.0 )

        # If PSF pixel size is different to map pixel size, rescale PSF accordingly
        if (pod['pix_arcsec']/psf_cdelt_arcsec)>1.001 or (pod['pix_arcsec']/psf_cdelt_arcsec)<0.999:
            zoom_factor = float(psf_cdelt_arcsec) / float(pod['pix_arcsec'])
            psf = scipy.ndimage.zoom(psf_in, (zoom_factor,zoom_factor), mode='nearest')

            # If necessary, rebin the PSF so it has odd dimensions
            if not ( psf.shape[0]%2!=0 and psf.shape[1]%2!=0 ):
                psf = ChrisFuncs.FromGitHub.martynbristow.rebin(psf, tuple(np.array(psf.shape)-1))

        # Else, if pixel sizes are already the same, leave as-is
        elif ((pod['pix_arcsec']/psf_cdelt_arcsec)>=0.999) and ((pod['pix_arcsec']/psf_cdelt_arcsec)<=0.001):
            psf = psf_in.copy()



    # Normalise PSF
    psf /= np.nansum(psf)

    # Background-subtract cutout
    cutout = pod['cutout'] - pod['bg_avg']

    # Produce mask for pixels we care about for fitting (ie, are inside photometric aperture and background annulus)
    mask = ChrisFuncs.Photom.EllipseMask(cutout, pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j']) #*band_dict['annulus_outer']

    # Produce guess values
    initial_sersic_amplitide = cutout[ int(round(pod['centre_i'])), int(round(pod['centre_j'])) ]
    initial_sersic_r_eff = pod['adj_semimaj_pix'] / 10.0
    initial_sersic_n = 1.0
    initial_sersic_x_0 = int(round(pod['centre_j']))
    initial_sersic_y_0 = int(round(pod['centre_i']))
    initial_sersic_ellip = ( pod['adj_axial_ratio'] - 1.0 ) / pod['adj_axial_ratio']
    initial_sersic_theta = np.deg2rad( pod['adj_angle'] )

    # Produce sersic model from guess parameters, for time trials
    sersic_x, sersic_y = np.meshgrid( np.arange(cutout.shape[1]), np.arange(cutout.shape[0]) )
    sersic_model = astropy.modeling.models.Sersic2D(amplitude=initial_sersic_amplitide, r_eff=initial_sersic_r_eff, n=initial_sersic_n, x_0=initial_sersic_x_0, y_0=initial_sersic_y_0, ellip=initial_sersic_ellip, theta=initial_sersic_theta)
    sersic_map = sersic_model(sersic_x,sersic_y)

    # Make sure that PSF array is smaller than sersic model array (as required for convolution); if not, remove its edges such that it is
    if psf.shape[0]>sersic_map.shape[0] or psf.shape[1]>sersic_map.shape[1]:
        excess = max( psf.shape[0]-sersic_map.shape[0], psf.shape[1]-sersic_map.shape[1] )
        border = max( 2, int( np.round( np.ceil( float(excess) / 2.0 ) - 1.0 ) ) )
        psf = psf[border:,border:]
        psf = psf[:-border,:-border]

    # Determine wither FFT convolution or direct convolution is faster for this kernel, using sersic model produced with guess parameters
    time_fft = time.time()
    conv_map = astropy.convolution.convolve_fft(sersic_map, psf, normalize_kernel=True)
    time_fft = time.time() - time_fft
    use_fft = True
    if not use_fft:
        if time_fft<10.0:
            use_fft = True
        else:
            time_direct = time.time()
            conv_map = astropy.convolution.convolve(sersic_map, psf, normalize_kernel=True)
            time_direct = time.time() - time_direct
            if time_fft<time_direct:
                use_fft = True
            else:
                use_fft = False



    # Set up parameters to fit galaxy with 2-dimensional sersic profile
    params = lmfit.Parameters()
    params.add('sersic_amplitide', value=initial_sersic_amplitide, vary=True)
    params.add('sersic_r_eff', value=initial_sersic_r_eff, vary=True, min=0.0, max=pod['adj_semimaj_pix'])
    params.add('sersic_n', value=initial_sersic_n, vary=True, min=0.1, max=10)
    params.add('sersic_x_0', value=initial_sersic_x_0, vary=False)
    params.add('sersic_y_0', value=initial_sersic_y_0, vary=False)
    params.add('sersic_ellip', value=initial_sersic_ellip, vary=True, min=0.5*initial_sersic_ellip, max=0.5*(1.0-initial_sersic_ellip)+initial_sersic_ellip)
    params.add('sersic_theta', value=initial_sersic_theta, vary=False)

    # Solve with LMfit to find parameters of best-fit sersic profile
    result = lmfit.minimize(Sersic_LMfit, params, args=(pod, cutout, psf, mask, use_fft), method='leastsq', ftol=1E-5, xtol=1E-5, maxfev=200)

    # Extract best-fit results
    sersic_amplitide = result.params['sersic_amplitide'].value
    sersic_r_eff = result.params['sersic_r_eff'].value
    sersic_n = result.params['sersic_n'].value
    sersic_x_0 = result.params['sersic_x_0'].value
    sersic_y_0 = result.params['sersic_y_0'].value
    sersic_ellip = result.params['sersic_ellip'].value
    sersic_theta = result.params['sersic_theta'].value

    # Construct underlying sersic map and convolved sersic map, using best-fit parameters
    sersic_model = astropy.modeling.models.Sersic2D(amplitude=sersic_amplitide, r_eff=sersic_r_eff, n=sersic_n, x_0=sersic_x_0, y_0=sersic_y_0, ellip=sersic_ellip, theta=sersic_theta)
    sersic_map = sersic_model(sersic_x, sersic_y)
    if use_fft==True:
        conv_map = astropy.convolution.convolve_fft(sersic_map, psf, normalize_kernel=True)
    elif use_fft==False:
        conv_map = astropy.convolution.convolve(sersic_map, psf, normalize_kernel=True)



    # Determine annulus properties before proceeding with photometry
    bg_inner_semimaj_pix = pod['adj_semimaj_pix'] * band_dict['annulus_inner']
    bg_width = (pod['adj_semimaj_pix'] * band_dict['annulus_outer']) - bg_inner_semimaj_pix

    # Evaluate pixels in source aperture and background annulus  unconvoled sersic map
    if float(band_dict['subpixel_factor'])==1.0:
        sersic_ap_calc = ChrisFuncs.Photom.EllipseSum(sersic_map, pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'])
        sersic_bg_calc = ChrisFuncs.Photom.AnnulusSum(sersic_map, bg_inner_semimaj_pix, bg_width, pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'])

    elif float(band_dict['subpixel_factor'])>1.0:
        sersic_ap_calc = ChrisFuncs.Photom.EllipseSumUpscale(sersic_map, pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'], upscale=band_dict['subpixel_factor'])
        sersic_bg_calc = ChrisFuncs.Photom.AnnulusSumUpscale(sersic_map, bg_inner_semimaj_pix, bg_width, pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'], upscale=band_dict['subpixel_factor'])

    # Background-subtract and measure unconvoled sersic source flux
    sersic_bg_clip = ChrisFuncs.SigmaClip(sersic_bg_calc[2], median=False, sigma_thresh=3.0)
    sersic_bg_avg = sersic_bg_clip[1] * float(band_dict['subpixel_factor'])**2.0
    #sersic_bg_avg = np.nanmedian(sersic_bg_calc[2])
    sersic_ap_sum = sersic_ap_calc[0] - (sersic_ap_calc[1] * sersic_bg_avg)

    # Evaluate pixels in source aperture and background annulus in convolved sersic map
    if float(band_dict['subpixel_factor'])==1.0:
        conv_ap_calc = ChrisFuncs.Photom.EllipseSum(conv_map, pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'])
        conv_bg_calc = ChrisFuncs.Photom.AnnulusSum(conv_map, bg_inner_semimaj_pix, bg_width, pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'])
    elif float(band_dict['subpixel_factor'])>1.0:
        conv_ap_calc = ChrisFuncs.Photom.EllipseSumUpscale(conv_map, pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'], upscale=band_dict['subpixel_factor'])
        conv_bg_calc = ChrisFuncs.Photom.AnnulusSumUpscale(conv_map, bg_inner_semimaj_pix, bg_width, pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'], upscale=band_dict['subpixel_factor'])

    # Background-subtract and measure convolved sersic source flux
    conv_bg_clip = ChrisFuncs.SigmaClip(conv_bg_calc[2], median=False, sigma_thresh=3.0)
    conv_bg_avg = conv_bg_clip[1] * float(band_dict['subpixel_factor'])**2.0
    #conv_bg_avg = np.nanmedian(conv_bg_calc[2])
    conv_ap_sum = conv_ap_calc[0] - (conv_ap_calc[1] * conv_bg_avg)



    # Find difference between flux measued on convoled and unconvoled sersic maps
    ap_correction = np.nanmax([ 1.0, (sersic_ap_sum/conv_ap_sum) ])

    # Apply aperture correction to pod, and return
    if kwargs_dict['verbose']: print '['+pod['id']+'] Applying aperture correction factor of '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(ap_correction,5))+'.'
    pod['ap_sum'] *= ap_correction
    pod['ap_error'] *= ap_correction
    return pod #astropy.io.fits.writeto('/home/saruman/spx7cjc/DustPedia/Conv.fits', conv_map, header=pod['in_header'], overwrite=True)





# Defie LMfit convolved-sersic function
def Sersic_LMfit(params, pod, cutout, psf, mask, use_fft, lmfit=True): #astropy.io.fits.writeto('/home/saruman/spx7cjc/DustPedia/Sersic.fits', sersic_map, overwrite=True)



    # Extract variable parameters
    sersic_amplitide = params['sersic_amplitide'].value
    sersic_r_eff = params['sersic_r_eff'].value
    sersic_n = params['sersic_n'].value
    sersic_x_0 = params['sersic_x_0'].value
    sersic_y_0 = params['sersic_y_0'].value
    sersic_ellip = params['sersic_ellip'].value
    sersic_theta = params['sersic_theta'].value

    # Generate sersic model given current paramters
    sersic_x, sersic_y = np.meshgrid( np.arange(cutout.shape[1]), np.arange(cutout.shape[0]) )
    sersic_model = astropy.modeling.models.Sersic2D(amplitude=sersic_amplitide, r_eff=sersic_r_eff, n=sersic_n, x_0=sersic_x_0, y_0=sersic_y_0, ellip=sersic_ellip, theta=sersic_theta)
    sersic_map = sersic_model( sersic_x, sersic_y )

    # Convolve sersic model with PSF
    if use_fft==True:
        conv_map = astropy.convolution.convolve_fft(sersic_map, psf, normalize_kernel=True)
    elif use_fft==False:
        conv_map = astropy.convolution.convolve(sersic_map, psf, normalize_kernel=True)

    # Calculate residuals, filtered by mask
    residuals = cutout - conv_map
    mask[ np.where(mask==0.0) ] = np.nan
    residuals *= mask
    residuals.flatten()
    residuals = residuals[ np.where( np.isnan(residuals)==False ) ]

    # Return residuals
    if lmfit==True:
        return residuals**2.0
    elif lmfit==False:
        return residuals, cutout-conv_map





# Define function that performs extinction correction on photometry, via ChrisFuncs function that calls IRSA dust extinction service (which uses the Schlafly & Finkbeiner 2011 prescription)
def ExtCorrrct(pod, source_dict, band_dict, kwargs_dict):



    # Run source details through function
    try:
        irsa_band_excorr = ChrisFuncs.ExtCorrrct(source_dict['ra'], source_dict['dec'], band_dict['band_name'], verbose=kwargs_dict['verbose'], verbose_prefix='['+pod['id']+'] ')
    except:
        if kwargs_dict['debug']==True:
            if kwargs_dict['verbose']: print '['+pod['id']+'] Extinction correction failed; entering debug mode.'
            pdb.set_trace()
        else:
            if kwargs_dict['verbose']: print '['+pod['id']+'] Extinction correction failed; skipping extinction correction.'
            irsa_band_excorr = 1

    # Update photometry with extinction corrections
    if 'ap_sum' in pod:  # photometry already done
        pod['ap_sum'] *= irsa_band_excorr
        pod['ap_error'] *= irsa_band_excorr
    elif kwargs_dict['save_images']:  # images as final product: modify cutout
        pod['cutout'] *= irsa_band_excorr

    # Report and return extinction-correced photometry
    irsa_band_excorr_mag = 2.51*np.log10(irsa_band_excorr)
    mag_str = 'NaN' if np.isnan(irsa_band_excorr_mag) else str(ChrisFuncs.FromGitHub.randlet.ToPrecision(irsa_band_excorr_mag,4))
    fac_str = 'NaN' if np.isnan(irsa_band_excorr) else str(ChrisFuncs.FromGitHub.randlet.ToPrecision(irsa_band_excorr,4))
    if np.isnan(irsa_band_excorr_mag): irsa_band_excorr_mag = 'NaN'
    if np.isnan(irsa_band_excorr): irsa_band_excorr = 'NaN'
    if kwargs_dict['verbose']: print '['+pod['id']+'] Applying Galactic extinction correction of '+mag_str+' mag (ie, factor of '+fac_str+').'
    return pod





# Define function that check is present band is to be excluded from photometry
def ExcludePhotom(pod, source_dict, band_dict, kwargs_dict):



    # Check if the photometry exclusion field actually contains characters; if so, make list of entries, and if not, record an empty list
    if isinstance(source_dict['photom_bands_exclude'], basestring):
        aperture_bands_exclude = source_dict['photom_bands_exclude'].split(';')
    elif source_dict['photom_bands_exclude']==False:
        aperture_bands_exclude = []

    # If present band is to be excluded, note this fact in pod
    if band_dict['band_name'] in aperture_bands_exclude:
        pod['band_exclude'] = True

    # If exclusion not required, record and return
    else:
        pod['band_exclude'] = False
        return pod

    # Create photometry output dictionry containing null values
    output_dict = {'band_name':band_dict['band_name'], 'ap_sum':np.NaN, 'ap_error':np.NaN}
    pod['null_output_dict'] = output_dict

    # Return pod
    if pod['verbose']: print '['+pod['id']+'] No photometry required from this source in this band.'

    return pod





# Define function to check that all bands which were supposed to undergo photometry have done so
def PhotomCheck(photom_attempts, photom_output_list, source_dict, bands_dict, kwargs_dict):



    # Compare number of bands with photometry returned to number of bands for which photometry was requested
    photom_limit = 3
    if len(photom_output_list)==len(bands_dict.keys()):
        photom_attempts = 'Complete'
        return photom_attempts, photom_output_list

    # Check how many attempts have been made so far, and proceed accordingly
    else:
        photom_attempts += 1
        time.sleep(30.0)
        if photom_attempts>=photom_limit:
            print '['+source_dict['name']+'] Photometry failed '+str(photom_limit)+' times in succession; suggest debugging.' #raise Exception('Photometry failed '+str(photom_limit)+' times in succession; suggest debugging.')
            if kwargs_dict['debug']:
                pdb.set_trace()

            # If debugging not enabled, record null photometry for missing bands
            else:
                photom_attempts = 'Complete'
                photom_output_list_revised = copy.deepcopy(photom_output_list)
                good_bands = []
                [ good_bands.append(photom_output['band_name']) for photom_output in photom_output_list ]
                all_bands = []
                for band in bands_dict.keys():
                    all_bands.append(bands_dict[band]['band_name'])
                for all_band in all_bands:
                    if all_band not in good_bands:
                        photom_output_list_revised.append( {'band_name':all_band, 'ap_sum':np.NaN, 'ap_error':np.NaN} )
                return photom_attempts, photom_output_list_revised

        # Increment up check and try again
        else:
            return photom_attempts, photom_output_list





# Define function that handles bands excluded from photometry, so that they appear in thumbnail grid
def ExcludedThumb(source_dict, bands_dict, kwargs_dict):



    # If thumbnails not required, end immediately
    if kwargs_dict['thumbnails']==False:
        return

    # Check if the photometry exclusion field for this source actually contains characters; if so make list of entries, else produce empty list
    if isinstance(source_dict['photom_bands_exclude'], basestring):
        photom_bands_exclude = source_dict['photom_bands_exclude'].split(';')
    else:
        photom_bands_exclude = []

    # Now consider bands which have been assigned a blancket photometry exclusion
    photom_bands_exclude = list( set( photom_bands_exclude ) )
    photom_bands_exclude = np.array(photom_bands_exclude)[ np.in1d( photom_bands_exclude, bands_dict.keys() ) ]

    # If no bands require processing here, end immediately; else prepare to loop over bands that do require processing
    if len(photom_bands_exclude)==0:
        return
    else:
        if kwargs_dict['verbose']: print '['+source_dict['name']+'] Preparing thumbnail data for bands excluded from photometry.'
        random.shuffle(photom_bands_exclude)

    # Find largest beam size and outer annulus size, for later usein thumbnail generation
    beam_arcsec_max = 0.0
    outer_annulus_max = 0.0
    pix_arcsec_max = 0.0
    for band_name in bands_dict:
        band_fitspath, band_file_found = CAAPR.CAAPR_Pipeline.FilePrelim(source_dict, bands_dict[band_name], kwargs_dict)
        if not band_file_found:#os.path.exists(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'.fits')):
            continue
        band_pix_matrix = astropy.wcs.WCS(astropy.io.fits.getheader(band_fitspath)).pixel_scale_matrix
        band_pix_arcsec = 3600.0 * np.sqrt( np.min(np.abs(band_pix_matrix))**2.0 + np.max(np.abs(band_pix_matrix))**2.0 )
        if band_pix_arcsec>pix_arcsec_max:
            pix_arcsec_max = band_pix_arcsec
        if bands_dict[band_name]['beam_arcsec']>beam_arcsec_max:
            beam_arcsec_max = bands_dict[band_name]['beam_arcsec']
        if bands_dict[band_name]['annulus_outer']>outer_annulus_max:
            outer_annulus_max = bands_dict[band_name]['annulus_outer']
    source_dict['beam_arcsec_max'] = beam_arcsec_max
    source_dict['outer_annulus_max'] = outer_annulus_max
    source_dict['pix_arcsec_max'] = pix_arcsec_max

    # In standard operation, process multiple sources in parallel
    if kwargs_dict['parallel']==True:
        ex_ap_pool = mp.Pool(processes=kwargs_dict['n_proc'])
        for band in photom_bands_exclude:
            ex_ap_pool.apply_async( ExcludedSubpipelinePhotom, args=(source_dict, bands_dict[band], kwargs_dict,) )
        ex_ap_pool.close()
        ex_ap_pool.join()
        del(ex_ap_pool)

    # If parallelisation is disabled, process sources one-at-a-time
    elif kwargs_dict['parallel']==False:
        for band in photom_bands_exclude:
            ExcludedSubpipelinePhotom(source_dict, bands_dict[band], kwargs_dict)





# Define 'pseudo-dummy' version of the photometry sub-pipeline, to run excluded bands through
def ExcludedSubpipelinePhotom(source_dict, band_dict, kwargs_dict_inviolate):
    source_id = source_dict['name']+'_'+band_dict['band_name']

    # Make deep copy of kwargs dict, to disable verbosity
    kwargs_dict = copy.deepcopy(kwargs_dict_inviolate)
    kwargs_dict['verbose'] = False

    # Run through initial stages of aperture sub-pipeline, as would occur usually
    in_fitspath_prelim, file_found = CAAPR.CAAPR_Pipeline.FilePrelim(source_dict, band_dict, kwargs_dict)
    if file_found == False:
        return
    pod = CAAPR.CAAPR_Pipeline.PodInitiate(in_fitspath_prelim, source_dict, band_dict, kwargs_dict)
    pod = CAAPR.CAAPR_Pipeline.MapPrelim(pod, source_dict, band_dict)
    if pod['within_bounds']==False:
        return
    pod = AperturePrelim(pod, source_dict, band_dict, kwargs_dict)
    pod = CAAPR.CAAPR_Pipeline.PolySub(pod, pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'], instant_quit=(not kwargs_dict['polysub']))
    CAAPR.CAAPR_IO.MemCheck(pod)

    # Use thumbnail cutout function to create a cutout that's only as large as it needs to be for the thumbnail grid
    semimaj_arcsec = np.sqrt( pod['adj_semimaj_arcsec']**2.0 - band_dict['beam_arcsec']**2.0 )
    thumb_rad_arcsec = np.ceil( 1.0 * source_dict['pix_arcsec_max'] ) + np.ceil( 1.75 * 0.5 * np.sqrt( (source_dict['outer_annulus_max']*2.0*semimaj_arcsec)**2.0 + (source_dict['beam_arcsec_max'])**2.0 ) )
    source_dict['thumb_rad_arcsec'] = thumb_rad_arcsec
    CAAPR.CAAPR_IO.ThumbCutout(source_dict, band_dict, kwargs_dict, pod['in_fitspath'], thumb_rad_arcsec)

    # Rename thumbnail cutout, and make it the 'active' map by repeating necessary processing
    thumb_output = os.path.join( kwargs_dict['temp_dir_path'], 'Processed_Maps', source_id+'_Thumbnail.fits' )
    pod['in_fitspath'] = thumb_output
    in_fitsdata = astropy.io.fits.open(pod['in_fitspath'])
    pod['in_image'] = in_fitsdata[0].data
    pod['in_header'] = in_fitsdata[0].header
    in_fitsdata.close()
    pod['in_wcs'] = astropy.wcs.WCS(pod['in_header'])
    pod['in_fitspath_size'] = float(os.stat(pod['in_fitspath']).st_size)
    thumb_centre_xy = pod['in_wcs'].wcs_world2pix( np.array([[ source_dict['ra'], source_dict['dec'] ]]), 0 )
    pod['centre_i'], pod['centre_j'] = float(thumb_centre_xy[0][1]), float(thumb_centre_xy[0][0])

    # Run thumbnail cutout thorugh AstroMagic (deleting any pre-existing data), save result, and remove temporary files
    if kwargs_dict['starsub']==True:
        pod['cutout'] = pod['in_image'].copy()
        pod['starsub_thumbnail'] = True
        if os.path.exists( os.path.join( kwargs_dict['temp_dir_path'], 'AstroMagic', band_dict['band_name'], source_dict['name']+'_'+band_dict['band_name']+'_StarSub.fits') ):
            os.remove( os.path.join( kwargs_dict['temp_dir_path'], 'AstroMagic', band_dict['band_name'], source_dict['name']+'_'+band_dict['band_name']+'_StarSub.fits') )
        pod = CAAPR.CAAPR_AstroMagic.Magic(pod, source_dict, band_dict, kwargs_dict)
        os.remove(thumb_output)
        magic_output = os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic', band_dict['band_name'], source_dict['name']+'_'+band_dict['band_name']+'_StarSub.fits')
        if os.path.exists(magic_output):
            os.remove(magic_output)
    else:
        pod['cutout'] = pod['in_image'].copy()
        os.remove(thumb_output)

    # Save resulting cutout
    astropy.io.fits.writeto(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_id+'.fits'), pod['cutout'], header=pod['in_header'], overwrite=True)
