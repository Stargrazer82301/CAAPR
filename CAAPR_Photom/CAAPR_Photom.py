# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import time
import math
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
import ChrisFuncs
import CAAPR_Pipeline
import CAAPR_IO
import CAAPR_AstroMagic
plt.ioff()





# The aperture-fitting sub-pipeline
def PipelinePhotom(source_dict, band_dict, kwargs_dict):
    source_id = source_dict['name']+'_'+band_dict['band_name']



    # Determine whether the user is specificing a directroy full of FITS files in this band (in which case use standardised filename format), or just a single FITS file
    try:
        os.path.isdir(band_dict['band_dir'])
    except:
        pdb.set_trace()
    if os.path.isdir(band_dict['band_dir']):
        in_fitspath = os.path.join( band_dict['band_dir'], source_dict['name']+'_'+band_dict['band_name'] )
    elif os.path.isfile(band_dict['band_dir']):
        in_fitspath = os.path.join( band_dict['band_dir'] )


    # Work out whether the file extension for FITS file in question is .fits or .fits.gz
    file_found = False
    if os.path.exists(in_fitspath+'.fits'):
        in_fitspath = in_fitspath+'.fits'
        file_found = True
    elif os.path.exists(in_fitspath+'.fits.gz'):
        in_fitspath = in_fitspath+'.fits.gz'
        file_found = True

    # Report to user if no file found; otherwise, progress to pipeline
    if file_found==False:
        print '['+source_id+'] No appropriately-named input file found in target directroy (please ensure that filesnames are in \"[NAME]_[BAND].fits\" format.)'
        print '['+source_id+'] Assuming no data in this band for current source.'
        output_dict = {'band_name':band_dict['band_name'],
                       'ap_sum':np.NaN,
                       'ap_error':np.NaN}
        return output_dict
    else:

        # Carry out small random wait, to stop RAM checks from syncing up
        time.sleep(5.0*np.random.rand())

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
               'band_dict':band_dict,
               'source_dict':source_dict,
               'kwargs_dict':kwargs_dict,
               'output_dir_path':kwargs_dict['output_dir_path'],
               'temp_dir_path':kwargs_dict['temp_dir_path'],
               'in_fitspath_size':in_fitspath_size,
               'id':source_id,
               'verbose':kwargs_dict['verbose']}
        CAAPR_IO.MemCheck(pod)



        # Read in aperture file
        aperture_table = np.genfromtxt(kwargs_dict['aperture_table_path'], delimiter=',', names=True, dtype=None)
        aperture_index = np.where( aperture_table['name']==source_dict['name'] )
        if aperture_index[0].shape[0]>1:
            raise ValueError('Aperture value caontains more than one entry for current galaxy')
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



        # Run pod through preliminary processing, to determine initial quantities; if target not within bounds of map, end processing here
        pod = CAAPR_Pipeline.MapPrelim(pod)
        if pod['within_bounds']==False:
            output_dict = {'band_name':band_dict['band_name'],
                       'ap_sum':np.NaN,
                       'ap_error':np.NaN}
            return output_dict



        # Adjust aperture to account for beam size of current band, and record to pod
        adj_semimaj_arcsec = 0.5 * ( (2.0*opt_semimaj_arcsec)**2.0 + band_dict['beam_arcsec']**2.0 )**0.5
        adj_semimin_arcsec = 0.5 * ( (2.0*opt_semimin_arcsec)**2.0 + band_dict['beam_arcsec']**2.0 )**0.5
        adj_axial_ratio = adj_semimaj_arcsec / adj_semimin_arcsec
        adj_angle = opt_angle
        pod['adj_semimaj_arcsec'] = adj_semimaj_arcsec
        pod['adj_semimin_arcsec'] = adj_semimin_arcsec
        pod['adj_semimaj_pix'] = adj_semimaj_arcsec / pod['pix_arcsec']
        pod['adj_semimin_pix'] = adj_semimin_arcsec / pod['pix_arcsec']
        pod['adj_axial_ratio'] = adj_axial_ratio
        pod['adj_angle'] = adj_angle



        # If star-removal is required, run pod through AstroMagic
        if kwargs_dict['starsub']==True:
            if band_dict['remove_stars']==True:
                if pod['verbose']: print '['+pod['id']+'] Removing foreground stars and background galaxies with PTS AstroMagic.'
                pod = CAAPR_AstroMagic.Magic(pod, source_dict, kwargs_dict, do_sat=False)



        # Run pod through function that removes large-scale sky using a 2-dimensional polynomial filter, with source aperture masked
        if kwargs_dict['polysub']==True:
            pod = CAAPR_Pipeline.PolySub(pod, pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'])



        # Run pod through function that (finally) performs the actual photometry
        pod = Photom(pod)



        # If photometry recorded null flux, skip determination of aperture noise, recording null value for this too, instead
        if np.isnan(pod['ap_sum'])==True:
            pod['ap_error'] = np.NaN
        else:



            # Attempt to determine aperture noise the preferred way, using full-size randomly-placed apertures
            if kwargs_dict['verbose']: print '['+source_id+'] Estimating aperture noise using full-size randomly-placed sky apertures.'
            ap_noise_dict = ApNoise(pod['cutout'].copy(), source_dict, band_dict, kwargs_dict, pod['adj_semimaj_pix'], pod['adj_axial_ratio'], pod['adj_angle'], pod['centre_i'], pod['centre_j'])
            if ap_noise_dict['fail']==False:
                if kwargs_dict['verbose']: print '['+source_id+'] Aperture noise successfully estimated using full-size randomly-placed sky apertures.'
                pod['ap_noise'] = ap_noise_dict['ap_noise']



            # If full-size apertures were unable to produce a valid aperture noise estimate, commence aperture extrapolation
            else:
                sky_success_counter = ap_noise_dict['sky_success_counter']
                if kwargs_dict['verbose']: print '['+source_id+'] Unable to estiamte aperture noise using full-size randomly-placed sky apertures (only '+str(int(sky_success_counter))+' could be placed); switching to aperture extrapolation.'
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
            if np.isnan(pod['ap_sum'])==False:
                if pod['verbose']:print '['+pod['id']+'] Final source flux and uncertainty is '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(pod['ap_sum'],4))+' +/- '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(pod['ap_error'],4))+' (in map units).'



        # If thumbnail images have been requested, save a copy of the current image (ie, with any star and/or background subtaction)
        if kwargs_dict['thumbnails']==True:
            astropy.io.fits.writeto(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_id+'.fits'), pod['cutout'], header=pod['in_header'])



        # Now return final photometry informaton to main pipeline, and clean up garbage
        output_dict = {'band_name':band_dict['band_name'],
                       'ap_sum':pod['ap_sum'],
                       'ap_error':pod['ap_error']}
        gc.collect()
        del(pod)
        return output_dict





# Define actual function that actually performs actual photometry
def Photom(pod):
    band_dict = pod['band_dict']
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

    # If more than 5% of the pixels inside the source aperture, or 75% of those in background annulus are NaN, record NaN flux, otherwise continue as per normal
    ap_nan_thresh = 0.05
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
def ApNoise(cutout, source_dict, band_dict, kwargs_dict, adj_semimaj_pix, adj_axial_ratio, adj_angle, centre_i, centre_j, mini=False):
    source_id = source_dict['name']+'_'+band_dict['band_name']
    cutout = cutout.copy()
    ap_debug = False



    # Handle input variables if mini-apertures are required
    if mini!=False:
        if isinstance(mini, float) or isinstance(mini, int):
            mini = float(mini)
            adj_semimaj_pix_full = adj_semimaj_pix
            adj_semimaj_pix = mini
        else:
            pdb.set_trace()
    else:
        adj_semimaj_pix_full = adj_semimaj_pix

    # Define charactaristics of circular aperture with same area as elliptical source aperture
    ap_area = np.pi * adj_semimaj_pix * ( adj_semimaj_pix / adj_axial_ratio )
    sky_ap_rad_pix = ( ap_area / np.pi )**0.5
    sky_border = int( sky_ap_rad_pix + 1.0 )#int( ( band_dict['annulus_outer'] * sky_ap_rad_pix ) + 1 )
    adj_semimin_pix = adj_semimaj_pix / adj_axial_ratio

    # Creating mask maps to describe no-go regions
    prior_mask = np.zeros(cutout.shape)
    #ChrisFuncs.Cutout(prior_mask, '/home/saruman/spx7cjc/DustPedia/Prior.fits')
    exclude_mask = ChrisFuncs.Photom.EllipseMask(cutout, adj_semimaj_pix_full, adj_axial_ratio, adj_angle, centre_i, centre_j)
    flag_mask = np.zeros(cutout.shape)
    attempt_mask = np.zeros(cutout.shape)

    # Set pixels in source aperture to all have NaN pixels, so they don't get sampled by sky annuli
    cutout_inviolate = cutout.copy()
    cutout[ np.where(ChrisFuncs.Photom.EllipseMask(cutout, adj_semimaj_pix_full, adj_axial_ratio, adj_angle, centre_i, centre_j)==1) ] = np.NaN



    # Define how many random aperture are desired/required/permitted
    sky_success_target = 50
    sky_success_min = 20
    sky_gen_max = 100

    # Generate random polar coordinates to draw from
    random_size = sky_success_target * sky_gen_max * 100
    random_failed = []
    random_theta_list = 360.0 * np.random.rand(random_size)
    random_r_list = (2.0*adj_semimin_pix) + np.abs(np.random.normal(loc=0.0, scale=5.0*adj_semimaj_pix_full, size=random_size))

    # Remove random coordinates that are more distant than most distant part of the coverage region the target source lies in
    cont_binary = np.zeros(cutout.shape)
    cont_binary[ np.where( np.isnan(cutout_inviolate)==False ) ] = 1
    cont_structure = np.array([[1,1,1], [1,1,1], [1,1,1]])
    cont_label = scipy.ndimage.measurements.label(cont_binary, structure=cont_structure)[0]
    cont_search_mask = ChrisFuncs.EllipseMask(cont_label, 3.0, 1.0, 0.0, centre_i, centre_j)
    cont_search_values = cont_label[ np.where( cont_search_mask==1 ) ]
    cont_target = scipy.stats.mode(cont_search_values[np.where(cont_search_values>0)])[0][0]
    cont_where_bad = np.where( cont_label!=cont_target )
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
    random_not_source = np.where( (abs(centre_i-random_i_list)>adj_semimaj_pix_full) & (abs(centre_j-random_j_list)>adj_semimaj_pix_full) )
    random_i_list = random_i_list[random_not_source]
    random_j_list = random_j_list[random_not_source]

    # Remove random coordinates that correspond to NaN pixels
    random_i_list_pix = np.round(random_i_list).astype(int)
    random_j_list_pix = np.round(random_j_list).astype(int)
    random_ij_values = cutout[(random_i_list_pix,random_j_list_pix)]
    random_ij_pix_good = np.where(np.isnan(random_ij_values)==False)
    random_i_list = random_i_list[random_ij_pix_good]
    random_j_list = random_j_list[random_ij_pix_good]





    # Commence creation of random sky apertures    sky_success_counter = 0
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
                if ap_debug: print 'Unable to generate suitable random sky aperture after '+str(sky_gen_max)+' attempts.'
                break



            # Select random coordinate set for this iteration
            random_accept = False
            while random_accept==False:
                random_index = int( np.round( np.random.rand() * float(random_i_list.shape[0]) ) )
                if random_index not in random_failed:
                    random_accept = True
            random_i = random_i_list[random_index]
            random_j = random_j_list[random_index]
            if ap_debug:
                ap_mask = ChrisFuncs.Photom.EllipseMask(cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
                attempt_mask[ np.where(ap_mask==1) ] = sky_success_counter
                print 'Aperture :'+str(sky_success_counter+1)+';   Generation: '+str(sky_gen_counter)+';   Pix Coords: ['+str(random_i)+','+str(random_j)+']'

            """
            # Produce r & theta coords for random sky aperture, dictated by gaussian distribution centred on source
            random_theta = 360.0 * np.random.rand()
            random_r = (2.0*adj_semimin_pix) + np.abs(np.random.normal(loc=0.0, scale=4.0*adj_semimaj_pix_full))
            random_i = centre_i + ( random_r * np.cos(np.radians(random_theta)) )#np.random.normal(loc=centre_i, scale=2.0*sky_ap_rad_pix)
            random_j = centre_j + ( random_r * np.sin(np.radians(random_theta)) )#np.random.normal(loc=centre_j, scale=2.0*sky_ap_rad_pix)

            # Check that generated sky aperture does not intersect border region
            if random_i<=sky_border or random_i>(cutout.shape[0]-sky_border) or random_j<=sky_border or random_j>(cutout.shape[1]-sky_border):
                if ap_debug: print 'Rejection: aperture in border'
                continue

            # Do crude basic that generated sky aperture does not intersect source; if it does, reject
            if ( abs(centre_i - random_i) < adj_semimaj_pix_full ) and ( abs(centre_j - random_j) < adj_semimaj_pix_full ):
                if ap_debug: print 'Rejection: aperture intersects source (according to basic check)'
                continue
            """
            # Do sophisticated check that generated sky aperture does not intersect source; if it does, reject
            exclude_sum = ChrisFuncs.Photom.EllipseSum(exclude_mask, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)[0]
            if exclude_sum>0:
                if ap_debug: print 'Rejection: aperture intersects source (according to sophisticated check)'
                random_failed.append(random_index)
                continue

            # Do basic chrck that the majority of the pixels in the generated sky aperture have not already been sampled by previous sky apertures; they have, reject
            prior_calc = ChrisFuncs.Photom.EllipseSum(prior_mask, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
            prior_calc[2][ np.where(prior_calc[2]>=1.0) ] = 1.0
            prior_frac = np.sum(prior_calc[2]) / float(prior_calc[1])
            if prior_frac>0.5:
                if ap_debug: print 'Rejection: aperture over-sampled (according to basic check)'
                random_failed.append(random_index)
                continue

            # Do sophisticated check that the majority of the pixels in the generated sky aperture have not already been sampled by previous sky apertures; they have, reject
            ap_mask_check = ChrisFuncs.Photom.EllipseMask(cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
            flag_mask_check = flag_mask.copy()
            flag_mask_check[np.where(ap_mask_check==1)] = int(2.0**(sky_success_counter+1.0))
            flag_tallies = np.array([ np.where(flag_mask_check==flag)[0].shape[0] for flag in (2.0**np.arange(0.0,sky_success_counter+2.0)).tolist() ])
            flag_check = np.where(flag_tallies<(0.5*ap_area))[0].shape[0]
            if flag_check>1:
                if ap_debug: print 'Rejection: aperture over-sampled (according to sophisticated check)'
                random_failed.append(random_index)
                continue

            # Evaluate pixels in sky aperture
            ap_calc = ChrisFuncs.Photom.EllipseSum(cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)

            # Evaluate pixels in sky annulus
            bg_inner_semimaj_pix = adj_semimaj_pix * band_dict['annulus_inner']
            bg_width = (adj_semimaj_pix * band_dict['annulus_outer']) - bg_inner_semimaj_pix
            bg_calc = ChrisFuncs.Photom.AnnulusSum(cutout, bg_inner_semimaj_pix, bg_width, 1.0, 0.0, random_i, random_j)

            # Check if more than 5% of the pixels inside the sky aperture, or 75% of those in sky annulus, are NaN; if so, reject
            try:
                if ap_calc[3][0].shape[0]==0 or ap_calc[1]==0:
                    ap_nan_frac = 0.0
                elif ap_calc[1]==0:
                    ap_nan_frac = 1.0
                else:
                    ap_nan_frac = float(ap_calc[3][0].shape[0]) / float(ap_calc[1]+float(ap_calc[3][0].shape[0]))
                if bg_calc[3][0].shape[0]==0:
                    bg_nan_frac = 0.0
                elif bg_calc[1]==0:
                    bg_nan_frac = 1.0
                else:
                    bg_nan_frac = float(bg_calc[3][0].shape[0]) / float(bg_calc[1]+bg_calc[3][0].shape[0])
            except:
                pdb.set_trace()
            #if ap_debug: print 'Aperture NaN fraction: '+str(ap_nan_frac)
            #if ap_debug: print 'Background NaN fraction: '+str(bg_nan_frac)
            ap_nan_thresh = 0.05
            bg_nan_thresh = 0.80
            if ap_nan_frac>ap_nan_thresh:
                if ap_debug: print 'Rejection: aperture contains too many NaNs'
                random_failed.append(random_index)
                continue
            if bg_nan_frac>bg_nan_thresh:
                if ap_debug: print 'Rejection: annulus contains too many NaNs'
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
                if ap_debug: print 'However, sufficient number of successful random apertures ('+str(int(sky_success_counter))+') already generated; proceeding'
                break

        # Calculate actual flux in sky aperture, and record
        bg_clip = ChrisFuncs.SigmaClip(bg_calc[2], median=False, sigma_thresh=3.0)
        bg_avg = bg_clip[1]
        ap_sum = ap_calc[0] - (ap_calc[1] * bg_avg)
        sky_sum_list.append(ap_sum)

        # Add this aperture to the prior mask and flag mask
        if not ap_debug: ap_mask = ChrisFuncs.Photom.EllipseMask(cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
        prior_mask += ap_mask
        flag_mask[np.where(ap_mask==1)] += 2.0**sky_success_counter

        # If target number of sky apertures have been processed, break out of loop
        if sky_success_counter>=sky_success_target:
            if ap_debug: print 'Target number of successful random apertures ('+str(int(sky_success_target))+') generated; proceeding'
            break



    # If total failure was encountered, end process and report now
    if sky_total_fail:
        ap_noise_dict = {'fail':True, 'prior_mask':prior_mask, 'flag_mask':flag_mask, 'sky_success_counter':sky_success_counter}
        cutout = cutout_inviolate
        return ap_noise_dict

    # Otherwise, calculate aperture noise using returned aperture values, and return
    else:
        sky_sum_list = np.array(sky_sum_list)
        ap_noise = ChrisFuncs.SigmaClip(sky_sum_list, tolerance=0.001, median=True, sigma_thresh=3.0)[0] #np.std(sky_sum_list)
        ap_noise_dict = {'fail':False, 'ap_noise':ap_noise, 'ap_num':sky_success_counter, 'prior_mask':prior_mask, 'flag_mask':flag_mask}
        #ChrisFuncs.Cutout(prior_mask, '/home/saruman/spx7cjc/DustPedia/Prior.fits')
        if kwargs_dict['verbose']:print '['+source_id+'] Aperture noise from current random apertures is '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(ap_noise,4))+' (in map units).'
        cutout = cutout_inviolate
        return ap_noise_dict





# Define funcion that attempts to estimate aperture noise using randomly-positioned sky apertures of given dimensions
def ApNoiseExtrap(cutout, source_dict, band_dict, kwargs_dict, adj_semimaj_pix, adj_axial_ratio, adj_angle, centre_i, centre_j):
    source_id = source_dict['name']+'_'+band_dict['band_name']



    # Define charactaristics of circular aperture with same area as elliptical source aperture
    ap_area = np.pi * adj_semimaj_pix * ( adj_semimaj_pix / adj_axial_ratio )
    sky_ap_rad_pix = ( ap_area / np.pi )**0.5

    # Generate list of mini-aperture sizes to use, and declare result lists
    mini_ap_rad_base = 2.0 #2.0**0.5 #1.25
    mini_ap_rad_pix_input = mini_ap_rad_base**np.arange( 1.0, np.ceil( 1.0 + math.log( sky_ap_rad_pix, mini_ap_rad_base ) ) )
    min_ap_rad_pix_output = []
    mini_ap_noise_output = []
    mini_ap_num_output = []

    # Loop over radii for mini-apertures
    for mini_ap_rad_pix in mini_ap_rad_pix_input:
        if kwargs_dict['verbose']:print '['+source_id+'] Finding aperture noise for mini-apertures of radius '+str(mini_ap_rad_pix)[:6]+' pixels.'
        mini_ap_noise_dict = ApNoise(cutout, source_dict, band_dict, kwargs_dict, adj_semimaj_pix, adj_axial_ratio, adj_angle, centre_i, centre_j, mini=mini_ap_rad_pix)

        # If mini-aperture succeeded, record and proceed; else, call it a day
        if mini_ap_noise_dict['fail']==False:
            min_ap_rad_pix_output.append(mini_ap_rad_pix)
            mini_ap_noise_output.append(mini_ap_noise_dict['ap_noise'])
            mini_ap_num_output.append(mini_ap_noise_dict['ap_num'])
        else:
            break

    # Convert output lists into arrays
    min_ap_rad_pix_output = np.array(min_ap_rad_pix_output)
    mini_ap_noise_output = np.array(mini_ap_noise_output)
    mini_ap_num_output = np.array(mini_ap_num_output).astype(float)
    #ChrisFuncs.Cutout(mini_ap_noise_dict['prior_mask'], '/home/saruman/spx7cjc/DustPedia/Prior.fits')



    # If insufficient points to make extrapolation, report failure; else proceed
    if min_ap_rad_pix_output.shape[0]<2:
        ap_noise_dict = {'fail':True, 'ap_noise':np.NaN}
        gc.collect()
        return ap_noise_dict
    else:



        # Calculate values to plot
        log_mini_ap_area = np.log10(np.pi*min_ap_rad_pix_output**2.0)
        log_mini_ap_noise = np.log10(mini_ap_noise_output)
        mini_ap_noise_err_rel = np.abs( log_mini_ap_noise * ( mini_ap_num_output**0.5 / mini_ap_num_output ) )
        mini_ap_noise_err = mini_ap_noise_err_rel * mini_ap_noise_output
        log_mini_ap_noise_err = ChrisFuncs.LogError(mini_ap_noise_output, mini_ap_noise_err)

        # Define straight-line function, and fit it to points
        def Line(x,m,c):
            return (m*x)+c
        line_fit = scipy.optimize.curve_fit(Line, log_mini_ap_area, log_mini_ap_noise, sigma=log_mini_ap_noise_err)
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





"""
 # Evaluate pixels in sky aperture; with consideration of sub-pixels if requested
ap_calc = ChrisFuncs.Photom.EllipseSum(cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
if float(band_dict['subpixel_factor'])==1.0:
    ap_calc = ChrisFuncs.Photom.EllipseSum(cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
elif float(band_dict['subpixel_factor'])>1.0:
    ap_calc = ChrisFuncs.Photom.EllipseSumUpscale(cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j, upscale=band_dict['subpixel_factor'])

# Evaluate pixels in sky annulus; with consideration of sub-pixels if requested
bg_inner_semimaj_pix = adj_semimaj_pix * band_dict['annulus_inner']
bg_width = (adj_semimaj_pix * band_dict['annulus_outer']) - bg_inner_semimaj_pix
bg_calc = ChrisFuncs.Photom.AnnulusSum(cutout, bg_inner_semimaj_pix, bg_width, 1.0, 0.0, random_i, random_j)
if float(band_dict['subpixel_factor'])==1.0:
    bg_calc = ChrisFuncs.Photom.AnnulusSum(cutout, bg_inner_semimaj_pix, bg_width, 1.0, 0.0, random_i, random_j)
elif float(band_dict['subpixel_factor'])>1.0:
    bg_calc = ChrisFuncs.Photom.AnnulusSumUpscale(cutout, bg_inner_semimaj_pix, bg_width, 1.0, 0.0, random_i, random_j, upscale=band_dict['subpixel_factor'])
"""