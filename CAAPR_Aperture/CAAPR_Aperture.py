# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import shutil
import pdb
import time
import copy
import numbers
import warnings
import psutil
import random
import multiprocessing as mp
import numpy as np
import scipy.optimize
import astropy.io.fits
import astropy.wcs
import astropy.convolution
import ChrisFuncs
import CAAPR





# The aperture-fitting sub-pipeline
def SubpipelineAperture(source_dict, band_dict, kwargs_dict):
    source_id = source_dict['name']+'_'+band_dict['band_name']



    # Carry out small random wait, to stop RAM checks from syncing up later
    time.sleep(5.0*np.random.rand())



    # Perform initial checks of target file type and location; return if not present
    in_fitspath, file_found = CAAPR.CAAPR_Pipeline.FilePrelim(source_dict, band_dict, kwargs_dict)
    if file_found == False:
        return



    # Create the pod (Photometry Organisation Dictionary), which will read in the FITS file, and bundle all the photometry data for this source & band into one dictionary to be passed between functions
    pod = CAAPR.CAAPR_Pipeline.PodInitiate(in_fitspath, source_dict, band_dict, kwargs_dict)



    # Run pod through preliminary processing, to determine initial quantities; if target not within bounds of map, end processing here
    pod = CAAPR.CAAPR_Pipeline.MapPrelim(pod, source_dict, band_dict)
    if pod['within_bounds']==False:
        return None
    CAAPR.CAAPR_IO.MemCheck(pod)



    # Check if this band is to be excluded from aperture-fitting; if so, return null aperture information
    pod = ExcludeAperture(pod, source_dict, band_dict, kwargs_dict)
    if pod['band_exclude']==True:
        return pod['null_output_dict']



    # If star-removal is required, run pod through AstroMagic
    pod = CAAPR.CAAPR_AstroMagic.Magic(pod, source_dict, band_dict, kwargs_dict)



    # Run pod through function that determines aperture shape, to provide preliminary estimate to facilitate removal of large-scale sky
    pod = ApertureShape(pod)



    # Run pod through function that removes large-scale sky using a 2-dimensional polynomial filter
    pod = CAAPR.CAAPR_Pipeline.PolySub( pod, 2.0*pod['semimaj_initial_pix'], pod['opt_axial_ratio'], pod['opt_angle'], instant_quit=max([not kwargs_dict['polysub'],pod['band_exclude']]) )



    # If sky polynomial removed, run pod through function that determines aperture shape, to provide final estiamte
    if pod['sky_poly']!=False:
        pod = ApertureShape(pod)



    # Run pod through function that determines aperture size
    pod = ApertureSize(pod, band_dict)



    # If thumbnail images have been requested, save a copy of the current image (ie, with any star and/or background subtaction)
    if kwargs_dict['thumbnails']==True:
        astropy.io.fits.writeto(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_id+'.fits'), pod['cutout'], header=pod['in_header'])



    # Now return final aperture informaton to main pipeline, and clean up garbage
    output_dict = {'band_name':band_dict['band_name'],
                   'opt_semimaj_arcsec':pod['opt_semimaj_arcsec'],
                   'opt_axial_ratio':pod['opt_axial_ratio'],
                   'opt_angle':pod['opt_angle']}
    gc.collect()
    del(pod)
    return output_dict





# Define function that determines the shape (not the size) of the source aperture in this band
def ApertureShape(pod):
    if pod['band_exclude']==True:
        return pod
    verbose = pod['verbose']
    if pod['verbose']: print '['+pod['id']+'] Commencing determination of appropriate axial ratio and positional angle for source aperture.'



    # Make preliminary per-pixel noise measurement by iteratively sigma-clipping cutout
    if verbose: print '['+pod['id']+'] Making preliminary per-pixel noise measurement.'
    clip_value = ChrisFuncs.SigmaClip(pod['cutout'], tolerance=0.001, sigma_thresh=3.0, median=True)
    noise_value = clip_value[0]
    field_value = clip_value[1]

    # Find all significant pixels that are connected to the region of the source (ie, withing a beam-width of the provided target coords)
    if verbose: print '['+pod['id']+'] Finding contiguous significant pixels around target.'
    semimaj_initial = int(round(pod['beam_pix']*1.0))
    cutoff = field_value + (4.0*noise_value)
    #ChrisFuncs.Cutout(cont_array_binary, '/home/saruman/spx7cjc/DustPedia/Cont.fits')
    #cont_structure = np.array([[1]*(2*semimaj_initial)]*(2*semimaj_initial))
    cont_array_prelim = ChrisFuncs.Photom.ContiguousPixels(pod['cutout'], semimaj_initial, pod['centre_i'], pod['centre_j'], cutoff)#, custom_structure=cont_structure)

    # Use binary erosion to remove thin artefacts (primarily diffraction spikes)
    erode_size = int(np.ceil(3.0*pod['beam_pix']))
    erode_centre = (0.5*float(erode_size))-0.5
    erode_structure = ChrisFuncs.Photom.EllipseMask(np.zeros([erode_size,erode_size]), pod['beam_pix'], 1.0, 0.0, erode_centre, erode_centre)
    erode_array = scipy.ndimage.morphology.binary_erosion(cont_array_prelim, structure=erode_structure).astype(int)
    #cont_array = scipy.ndimage.measurements.label(erode_array)[0]
    cont_array = ChrisFuncs.Photom.ContiguousPixels(erode_array, semimaj_initial, pod['centre_i'], pod['centre_j'], 1E-50)

    # If remainging contiguous pixel region has same or fewer number of pixels than erosion structure, replace with erosion sturcture
    if np.sum(cont_array)<=(np.sum(erode_structure)):
        cont_array = erode_structure

    # Find ellipse that best fits outline of contiguous region
    if verbose: print '['+pod['id']+'] Fitting ellipse to perimeter of contiguous significant pixels.'
    cont_x = ((np.where(cont_array==1))[1])
    cont_y = ((np.where(cont_array==1))[0])
    if cont_x.shape[0]>10:
        try:
            cont_ellipse = ChrisFuncs.Photom.EllipseFit(cont_x, cont_y)
            opt_axial_ratio = max(cont_ellipse[1]) / min(cont_ellipse[1])
            opt_angle = cont_ellipse[2]
            semimaj_initial = max(cont_ellipse[1])
        except:
            opt_axial_ratio = 1.0
            opt_angle = 0.0

    # If too few significant pixels, default to circular aperture
    else:
        opt_axial_ratio = 1.0
        opt_angle = 0.0
    if verbose: print '['+pod['id']+'] Ellipse angle: '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(opt_angle,4))+' degrees; Ellipse axial ratio: '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(opt_axial_ratio,4))+'.'

    # Clean garbage, record results to pod, and return
    gc.collect()
    pod['cutout_clip'] = clip_value
    pod['opt_axial_ratio'] = opt_axial_ratio
    pod['opt_angle'] = opt_angle
    pod['semimaj_initial_pix'] = semimaj_initial
    pod['semimaj_initial_arcsec'] = semimaj_initial * pod['pix_arcsec']
    return pod





# Define function that determines the size of the source aperture in this band
def ApertureSize(pod, band_dict):
    if pod['band_exclude']:
        return pod
    if pod['verbose']: print '['+pod['id']+'] Commencing determination of appropriate size for source aperture.'
    verbose = pod['verbose']

    # Define sub-function that determines SNR of a defined annulus; and if requested, determines residual between the SNR of a defined annulus, and a target SNR of 2
    def AnnulusSNR(semimaj, pod, cutout, width, i_trans, j_trans, residual):
        semimaj = semimaj[0]
        sig_annulus = ChrisFuncs.Photom.AnnulusQuickSum(cutout, semimaj, width, pod['opt_axial_ratio'], pod['opt_angle'], pod['centre_i'], pod['centre_j'], i_trans, j_trans)
        sig_value = ChrisFuncs.SigmaClip(sig_annulus[2], tolerance=0.005, sigma_thresh=2.0, median=True)[1]#np.median(sig_annulus[2])
        noise_value = pod['cutout_clip'][0]
        field_value = pod['cutout_clip'][1]
        ann_SNR = (sig_value - field_value) / noise_value
        if residual==True:
            ann_residual = abs(2.0-ann_SNR)
            #print 'SNR: '+str(ann_SNR)+', Semi-Maj: '+str(semimaj)+', Residual:'+str(ann_residual)
            return ann_residual
        elif residual==False:
            #print 'SNR: '+str(ann_SNR)+', Semi-Maj: '+str(semimaj)
            return ann_SNR



    # Construct kernel with FWHM equal to 3 beam-widths, by which to smooth map
    if verbose: print '['+pod['id']+'] Convolving map to lower resolution (twice the beam width) for radial analysis.'
    pix_size = pod['pix_arcsec']
    res_in = band_dict['beam_arcsec']
    res_out = 2.0*band_dict['beam_arcsec']#36.0
    kernel_fwhm = np.sqrt( (res_out/pix_size)**2.0 - (res_in/pix_size)**2.0 )

    # Determine if map contains NaN pixels, excluding those that simply represent edge of map
    cutout_prelabel = pod['cutout'].copy()
    cutout_prelabel = cutout_prelabel.byteswap().newbyteorder().astype('float64')
    cutout_prelabel[ np.where(np.isnan(cutout_prelabel)) ] = 0.0
    cutout_label = scipy.ndimage.label(cutout_prelabel)
    cutout_preconv = pod['cutout'].copy()
    cutout_preconv[ np.where(cutout_label==0) ] = 0.0

    # If map contains no NaN pixels, smooth it the quick Scipy way
    cutout_unconv = pod['cutout'].copy()
    if np.where(np.isnan(cutout_preconv)==True)[0].shape[0]==0:
        if verbose: print '['+pod['id']+'] No NaN pixels within coverage area; convolving using quick method.'
        pod['cutout'] = scipy.ndimage.filters.gaussian_filter(cutout_preconv, kernel_fwhm)

    # Else if map contains NaNs, do it the robust (but very-very slow, very-very memory intensive) Astropy way
    else:
        if verbose: print '['+pod['id']+'] NaN pixels within coverage area; convolving using slower NaN-compatible method.'
        CAAPR.CAAPR_IO.MemCheck(pod, thresh_factor=20.0)
        kernel = astropy.convolution.kernels.Gaussian2DKernel(kernel_fwhm)
        pod['cutout'] = astropy.convolution.convolve_fft(pod['cutout'], kernel, interpolate_nan=False, normalize_kernel=True, ignore_edge_zeros=False, allow_huge=True)
        pod['cutout'][ np.where( np.isnan(cutout_unconv)==True ) ] = np.NaN



    # Prepare arrays of transposed coordinates, to allow for rapid radial evaluating
    if verbose: print '['+pod['id']+'] Constructing arrays of transposed radial coordinates.'
    coords_trans = ChrisFuncs.Photom.AnnulusQuickPrepare(pod['cutout'], pod['opt_angle'], pod['centre_i'], pod['centre_j'])
    i_trans, j_trans = coords_trans[0], coords_trans[1]

    # To start with, to make new estimate of map noise that isn't contaminated by the target galaxy, by masking all pixels beyond semi-major axis suggested by contiguous significant pixels.
    brute_mask = ChrisFuncs.Photom.EllipseMask(pod['cutout'], pod['semimaj_initial_pix'], pod['opt_axial_ratio'], pod['opt_angle'], pod['centre_i'], pod['centre_j'])
    cutout_brute_masked = pod['cutout'].copy()
    cutout_brute_masked[ np.where( brute_mask==1 ) ] = np.nan
    cutout_clip_masked = ChrisFuncs.SigmaClip(cutout_brute_masked, tolerance=0.001, sigma_thresh=3.0, median=True)
    pod['cutout_clip'] = cutout_clip_masked

    # Now, perform a coarse brute force ckeck of a small number of radii over a wide range, to find rough location of edge of the source
    if verbose: print '['+pod['id']+'] Finding size of target source with coarse analysis.'
    ann_brute_range = np.linspace(0.75*pod['semimaj_initial_pix'], 2.25*pod['semimaj_initial_pix'], num=15)
    ann_brute_range = ann_brute_range[::-1]
    ann_brute_width = abs( ann_brute_range[1] - ann_brute_range[0] )
    snr_success = False
    for i in range(0, len(ann_brute_range)):
        ann_brute_snr = AnnulusSNR([ann_brute_range[i]], pod, pod['cutout'], ann_brute_width, i_trans, j_trans, False)
        if ann_brute_snr>2:
            snr_success = True
            #ann_brute_semimaj = ann_brute_range[i-1]
            ann_bounds = [( ann_brute_range[ max(i-2,0) ], ann_brute_range[ min(i+1,len(ann_brute_range)-1) ] )]
            if verbose: print '['+pod['id']+'] Course analysis finds that radial SNR=2 between semi-major axes of '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(ann_bounds[0][1]*pod['pix_arcsec'],4))+' and '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(ann_bounds[0][0]*pod['pix_arcsec'],4))+' arcseconds.'
            break

    # If SNR=2 threshold not reached, set to minimum semi-major axis of two beam-widths
    if snr_success==False:
        ann_bounds = [( ann_brute_range[i], np.floor(pod['cutout'].shape[0]/2.0) )]
        opt_semimaj_pix = pod['beam_pix'] * 2.0
        opt_semimaj_arcsec = opt_semimaj_pix * pod['pix_arcsec']
        if verbose: print '['+pod['id']+'] No SNR=2 threshold found; hence reverting to two beam-width minimum value of '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(opt_semimaj_arcsec,4))+' arcseconds.'
    else:

        # Now use scipy differential evolution optimisation to find, with more precision, the semi-major axis at which annulus falls to a SNR of 2
        ann_beams = 1.0
        ann_width = np.ceil( ann_beams * pod['beam_pix'] )
        if verbose: print '['+pod['id']+'] Refining size of target source with more precise analysis.'
        ann_fit = scipy.optimize.differential_evolution(AnnulusSNR, ann_bounds, args=(pod, pod['cutout'], ann_width, i_trans, j_trans, True), maxiter=5, popsize=10, polish=False)
        """
        ann_guess = ann_brute_semimaj
        ann_fit = scipy.optimize.minimize(AnnulusSNR, ann_guess, args=(pod, pod['cutout'], ann_width, i_trans, j_trans))#method='Nelder-Mead', tol=1E-4, method='L-BFGS-B', bounds=[(pod['semimaj_initial_pix'], None)],
        minimizer_kwargs = {'args':(pod, pod['cutout'], ann_width, i_trans, j_trans)}
        ann_fit = scipy.optimize.basinhopping(AnnulusSNR, ann_guess, T=0.1, stepsize=5.0, minimizer_kwargs=minimizer_kwargs)
        """
        # Extract results from fitting
        opt_semimaj_pix = ann_fit['x'][0]
        opt_semimaj_arcsec = opt_semimaj_pix * pod['pix_arcsec']
        if verbose: print '['+pod['id']+'] Precision analysis finds that radial SNR=2 at semi-major axis of '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(opt_semimaj_arcsec,4))+' arcseconds.'

        # For small sources, default to minimum semi-major axis of two beam-widths
        if opt_semimaj_pix<(ann_beams*2.0):
            opt_semimaj_pix = pod['beam_pix'] * 2.0
            opt_semimaj_arcsec = opt_semimaj_pix * pod['pix_arcsec']
            if verbose: print '['+pod['id']+'] Semi-major axis at which SNR=2 is less than two beam-widths; hence reverting to two beam-width minimum value of '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(opt_semimaj_arcsec,4))+' arcseconds.'

    # Establish what fraction of the pixels inside a band's aperture are NaNs
    pix_good = ChrisFuncs.Photom.EllipseQuickSum(pod['cutout'], opt_semimaj_pix, pod['opt_axial_ratio'], pod['opt_angle'], pod['centre_i'], pod['centre_j'], i_trans, j_trans)[1]
    pix_tot = np.where( ChrisFuncs.Photom.EllipseMask(pod['cutout'], opt_semimaj_pix, pod['opt_axial_ratio'], pod['opt_angle'], pod['centre_i'], pod['centre_j']) == 1 )[0].shape[0]
    pdb.set_trace()
    if pix_tot==0.0:
        pix_good_frac = 0.0
    else:
        pix_good_frac = float(pix_good) / float(pix_tot)

    # Before final reporting, tidy up and return to unconvolved cutout
    pod['cutout'] = cutout_unconv
    del(cutout_unconv)
    gc.collect()

    # If more than 10% of the pixels in the aperture are NaNs, report NaN values for aperture dimensions; else proceed normally
    if pix_good_frac<0.9:
        opt_semimaj_pix = pod['beam_pix'] * 2.0
        opt_semimaj_arcsec = opt_semimaj_pix * pod['pix_arcsec']
        if verbose: print '['+pod['id']+'] More than 10% of pixels in fitted aperture are NaNs; hence reverting to two beam-width minimum value of '+str(opt_semimaj_arcsec)[:7]+' arcseconds.'
        pod['opt_semimaj_arcsec'] = opt_semimaj_arcsec
        pod['opt_axial_ratio'] = 1.0
        pod['opt_angle'] = 0.0
    else:

        # Deconvolve aperture semi-major axis with beam, by subtracting in quadrature
        adj_semimaj_arcsec = abs( opt_semimaj_arcsec**2.0 - (0.5*band_dict['beam_arcsec'])**2.0 )**0.5
        opt_semimin_arcsec = opt_semimaj_arcsec / pod['opt_axial_ratio']
        adj_semimin_arcsec = abs( opt_semimin_arcsec**2.0 - (0.5*band_dict['beam_arcsec'])**2.0 )**0.5
        adj_ax_ratio = adj_semimaj_arcsec / adj_semimin_arcsec

        # Record final dimensions to pod
        pod['opt_semimaj_arcsec'] = adj_semimaj_arcsec
        pod['opt_semimaj_pix'] = adj_semimaj_arcsec / pod['pix_arcsec']
        pod['opt_axial_ratio'] = adj_ax_ratio

    # Clean up, then return results
    gc.collect()
    return pod





# Define function that combines a set of apertures for a given source into
def CombineAperture(aperture_output_list, source_dict, kwargs_dict):
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Combining individual apertures from all bands to generate final aperture.'



    # Extract various aperture values
    semimaj_arcsec_list = []
    axial_ratio_list = []
    angle_list = []
    for aperture in aperture_output_list:
        try:
            semimaj_arcsec_list.append( aperture['opt_semimaj_arcsec'] )
            axial_ratio_list.append( aperture['opt_axial_ratio'] )
            angle_list.append( aperture['opt_angle'] )
        except:
            pdb.set_trace()

    # Find largest semi-major axis, and use to define size of enclosisity array (which will have pixels some fraction the size of the smallest semi-major axis)
    semimaj_max = np.nanmax(semimaj_arcsec_list) #semimaj_min = np.nanmin(semimaj_arcsec_list)
    ap_array_pix_size = 0.005 * semimaj_max
    ap_array_scale = int( np.round( semimaj_max / ap_array_pix_size ) )
    ap_array = np.zeros([ 1+(2.2*ap_array_scale), 1+(2.2*ap_array_scale) ])
    centre_i, centre_j = 1+(1.1*ap_array_scale), 1+(1.1*ap_array_scale)
    semimaj_pix_list = np.array(semimaj_arcsec_list) / ap_array_pix_size

    # Loop over each aperture, adding to enclosisity array
    for a in range(0, len(semimaj_pix_list)):
        if np.isnan(semimaj_pix_list[a])==False:
            ap_mask = ChrisFuncs.EllipseMask(ap_array, semimaj_pix_list[a], axial_ratio_list[a], angle_list[a], centre_i, centre_j)
            ap_array[ np.where( ap_mask==1 ) ] += 1
    #ChrisFuncs.Cutout(ap_array, '/home/saruman/spx7cjc/DustPedia/Ap.fits')

    # Find ellipse that traces edge of enclosisity region
    cont_rad_initial_pix = 2.0#( semimaj_min / np.nanmax(axial_ratio_list) ) / ap_array_pix_size
    cont_array = ChrisFuncs.Photom.ContiguousPixels(ap_array, cont_rad_initial_pix, centre_i, centre_j, 0.1)
    cont_x = ((np.where(cont_array==1))[1])
    cont_y = ((np.where(cont_array==1))[0])
    if cont_x.shape[0]>10:
        try:
            cont_ellipse = ChrisFuncs.Photom.EllipseFit(cont_x, cont_y)
            cont_axial_ratio = max(cont_ellipse[1]) / min(cont_ellipse[1])
            cont_angle = cont_ellipse[2]
            cont_semimaj_pix = max([ cont_ellipse[1].max(), np.nanmax(semimaj_pix_list) ])
        except:
            cont_axial_ratio = 1.0
            cont_angle = 0.0
            cont_semimaj_pix = max([ cont_ellipse[1].max(), np.nanmax(semimaj_pix_list) ])
    else:
        pdb.set_trace()
        cont_axial_ratio = 1.0
        cont_angle = 0.0
        cont_semimaj_pix = 2.0 * np.max('beam_width')

    # Convert final semi-major axis back to arcsec and apply expanson factor, then clean garbage and return results
    if isinstance(kwargs_dict['expansion_factor'], float) or isinstance(kwargs_dict['expansion_factor'], int):
        expansion_facor = float(kwargs_dict['expansion_factor'])
    else:
        expansion_facor = 1.0
    cont_semimaj_arcsec = cont_semimaj_pix * ap_array_pix_size * expansion_facor

    # If final aperture is smaller than defined minimum aperture, switch to defined minimum
    if cont_semimaj_arcsec<source_dict['fitting_min_semimaj_arcsec']:
        if kwargs_dict['verbose']: print '['+source_dict['name']+'] Fitted aperture is smaller than minimum permitted aperture size; reverting to minimum permitted aperture size.'
        cont_semimaj_arcsec = source_dict['fitting_min_semimaj_arcsec']
        cont_axial_ratio = 1.0
        cont_angle = 0.0

    # Clean garbage and return results
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Final ellipse semi-major axis: '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(cont_semimaj_arcsec,4))+' arcsec; final ellipse angle: '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(cont_angle,4))+' '+' degrees; final ellipse axial ratio: '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(cont_axial_ratio,4))+'.'
    gc.collect()
    return [cont_semimaj_arcsec, cont_axial_ratio, cont_angle, ap_array]





# Define function that check is present band is to be excluded from aperture-fitting
def ExcludeAperture(pod, source_dict, band_dict, kwargs_dict):



    # Check if the aperture exclusion field actually contains characters; if so, make list of entries, and if not, record an empty list
    if isinstance(source_dict['aperture_bands_exclude'], str):
        aperture_bands_exclude = source_dict['aperture_bands_exclude'].split(';')
    elif source_dict['aperture_bands_exclude']==False:
        aperture_bands_exclude = []

    # If present band is to be excluded, note this fact in pod
    if (band_dict['consider_aperture']==False) or (band_dict['band_name'] in aperture_bands_exclude):
        pod['band_exclude'] = True

    # If exclusion not required, record and return
    else:
        pod['band_exclude'] = False
        return pod

    # Set generic null aperture properties
    pod['opt_axial_ratio'] = 1.0
    pod['opt_angle'] = 0.0
    pod['opt_semimaj_arcsec'] = ( (2.0*band_dict['beam_arcsec'])**2.0 - band_dict['beam_arcsec']**2.0 )**0.5
    pod['opt_semimaj_pix'] = pod['opt_semimaj_arcsec'] / pod['pix_arcsec']
    pod['semimaj_initial_pix'] = pod['opt_semimaj_arcsec'] / pod['pix_arcsec']

    # Create aperture output dictionry containing null values
    output_dict = {'band_name':band_dict['band_name'],
                   'opt_semimaj_arcsec':pod['opt_semimaj_arcsec'],
                   'opt_axial_ratio':pod['opt_axial_ratio'],
                   'opt_angle':pod['opt_angle']}
    pod['null_output_dict'] = output_dict

    # Return pod
    if pod['verbose']: print '['+pod['id']+'] No aperture fitting required from this source in this band.'
    return pod





# Define function that handles bands excluded from aperture fitting, so that they appear in thumbnail grid
def ExcludedThumb(source_dict, bands_dict, kwargs_dict, aperture_list, aperture_combined):



    # If thumbnails not required, end immediately
    if kwargs_dict['thumbnails']==False:
        return

    # Check if the aperture exclusion field for this source actually contains characters; if so make list of entries, else produce empty list
    if isinstance(source_dict['aperture_bands_exclude'], str):
        aperture_bands_exclude = source_dict['aperture_bands_exclude'].split(';')
    else:
        aperture_bands_exclude = []

    # Now consider bands which have been assigned a blancket aperture exclusion
    [ aperture_bands_exclude.append(band) for band in bands_dict.keys() if bands_dict[band]['consider_aperture']==False ]
    aperture_bands_exclude = list( set( aperture_bands_exclude ) )
    aperture_bands_exclude = np.array(aperture_bands_exclude)[ np.in1d( aperture_bands_exclude, bands_dict.keys() ) ]

    # If no bands require processing here, end immediately; else prepare to loop over bands that do require processing
    if len(aperture_bands_exclude)==0:
        return
    else:
        if kwargs_dict['verbose']: print '['+source_dict['name']+'] Preparing thumbnail data for bands excluded from aperture-fitting.'
        random.shuffle(aperture_bands_exclude)

    # In standard operation, process multiple sources in parallel
    if kwargs_dict['parallel']==True:
        ex_ap_pool = mp.Pool(processes=kwargs_dict['n_proc'])
        for band in aperture_bands_exclude:
            ex_ap_pool.apply_async( ExcludedSubpipelineAperture, args=(aperture_combined, source_dict, bands_dict[band], kwargs_dict,) )
        ex_ap_pool.close()
        ex_ap_pool.join()
        del(ex_ap_pool)

    # If parallelisation is disabled, process sources one-at-a-time
    elif kwargs_dict['parallel']==False:
        for band in aperture_bands_exclude:
            ExcludedSubpipelineAperture(aperture_combined, source_dict, bands_dict[band], kwargs_dict)





# Define 'pseudo-dummy' version of the aperture sub-pipeline, to run excluded bands through
def ExcludedSubpipelineAperture(aperture_combined, source_dict, band_dict, kwargs_dict_inviolate):
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
    CAAPR.CAAPR_IO.MemCheck(pod)

    # Use thumbnail cutout function to create a cutout that's only as large as it needs to be for the thumbnail grid
    img_naxis_pix = np.max([ pod['in_header']['NAXIS1'], pod['in_header']['NAXIS1'] ])
    img_naxis_arcsec = float(img_naxis_pix) * float(pod['pix_arcsec'])
    img_rad_arcsec = img_naxis_arcsec / 2.0
    thumb_rad_arcsec = 1.25 * aperture_combined[0] * band_dict['annulus_outer']
    CAAPR.CAAPR_IO.ThumbCutout(source_dict, band_dict, kwargs_dict, pod['in_fitspath'], img_rad_arcsec, thumb_rad_arcsec)

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

    # Run thumbnail cutout thorugh AstroMagic, save result, and delete temporary files
    pod['cutout'] = pod['in_image'].copy()
    pod['starsub_thumbnail'] = True
    pod = CAAPR.CAAPR_AstroMagic.Magic(pod, source_dict, band_dict, kwargs_dict)
    os.remove(thumb_output)
    astropy.io.fits.writeto(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_id+'.fits'), pod['cutout'], header=pod['in_header'], clobber=True)
    magic_output = os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic', band_dict['band_name'], source_dict['name']+'_'+band_dict['band_name']+'_StarSub.fits')
    if os.path.exists(magic_output):
        os.remove(magic_output)






