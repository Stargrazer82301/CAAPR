# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import time
import numbers
import warnings
import psutil
import numpy as np
import scipy.optimize
import astropy.io.fits
import astropy.wcs
import astropy.convolution
import ChrisFuncs
import CAAPR_Pipeline
import CAAPR_IO





# The aperture-fitting sub-pipeline
def PipelineAperture(source_dict, band_dict, output_dir_path, temp_dir_path, verbose):
    source_id = source_dict['name']+'_'+band_dict['band_name']



    # Using standard filename format, construct full file path, and work out whether the file extension is .fits or .fits.gz
    in_fitspath = os.path.join( band_dict['band_dir'], source_dict['name']+'_'+band_dict['band_name'] )
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
        #raise ValueError('No appropriately-named input file found in target directroy (please ensure that filesnames are in \"[NAME]_[BAND].fits\" format.')
    else:

        # Before processing, check enough RAM is free
        in_fitspath_size = float(os.stat(in_fitspath).st_size)
        mem_wait = True
        while mem_wait:
            mem_wait = CAAPR_IO.MemCheck(in_fitspath_size, verbose=verbose)
            if mem_wait:
                print '['+source_id+'] Waiting for some RAM to free up before commencing processing.'
                time.sleep(10.0+(5.0*np.random.rand()))
            elif not mem_wait:
                break

        # Read in FITS file in question
        in_fitsdata = astropy.io.fits.open(in_fitspath)
        in_image = in_fitsdata[0].data
        in_header = in_fitsdata[0].header
        in_fitsdata.close()
        in_wcs = astropy.wcs.WCS(in_header)

        # Create the pod (Photometry Organisation Dictionary), which will bundle all the photometry data for this source & band into one dictionary to be passed between functions
        pod = {'in_fitspath':in_fitspath,
               'in_image':in_image,
               'in_header':in_header,
               'in_wcs':in_wcs,
               'cutout':in_image.copy(),
               'band_dict':band_dict,
               'source_dict':source_dict,
               'output_dir_path':output_dir_path,
               'temp_dir_path':temp_dir_path,
               'in_fitspath_size':in_fitspath_size,
               'id':source_id,
               'verbose':verbose}

        # Run pod through preliminary processing, to determine initial quantities; if target not within bounds of map, end processing here
        pod = CAAPR_Pipeline.MapPrelim(pod)
        if pod['within_bounds']==False:
            return pod

        # If star-removal is required, run pod through AstroMagic
        if band_dict['remove_stars']==True:
            CAAPR_Pipeline.AstroMagic(pod)

        # Run pod through function that determines aperture shape, to provide preliminary estimate to facilitate removal of large-scale sky
        pod = ApertureShape(pod)

        # Run pod through function that removes large-scale sky using a 2-dimensional polynomial filter
        pod = CAAPR_Pipeline.PolySub(pod)

        # If sky polynomial removed, run pod through function that determines aperture shape, to provide final estiamte
        if pod['sky_poly']!=False:
            pod = ApertureShape(pod)

        # Run pod through function that determines aperture size
        pod = ApertureSize(pod)

        # Now return final aperture informaton to main pipeline
        output_dict = {'band_name':band_dict['band_name'],
                       'opt_semimaj_arcsec':pod['opt_semimaj_arcsec'],
                       'opt_axial_ratio':pod['opt_axial_ratio'],
                       'opt_angle':pod['opt_angle']}
        return output_dict





# Define function that determines the shape (not the size) of the source aperture in this band
def ApertureShape(pod, verbose=False):
    if verbose: print '['+pod['id']+'] Commencing determination of appropriate axial ratio and positional angle for source aperture.'
    verbose = pod['verbose']



    # Make preliminary per-pixel noise measurement by iteratively sigma-clipping cutout
    if verbose: print '['+pod['id']+'] Making preliminary per-pixel noise measurement.'
    clip_value = ChrisFuncs.SigmaClip(pod['cutout'], tolerance=0.001, sigma_thresh=3.0, median=True)
    noise_value = clip_value[0]
    field_value = clip_value[1]

    # Find all significant pixels that are connected to the region of the source (ie, withing a beam-width of the provided target coords)
    if verbose: print '['+pod['id']+'] Finding contiguous significant pixels around target.'
    semimaj_initial = int(round(pod['beam_pix']*1.0))
    cutoff = field_value + (3.0*noise_value)
    #ChrisFuncs.Cutout(cont_array_binary, '/home/saruman/spx7cjc/DustPedia/Cont.fits')
    cont_array = ChrisFuncs.Photom.ContiguousPixels(pod['cutout'], semimaj_initial, pod['centre_i'], pod['centre_j'], cutoff)

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
    if verbose: print '['+pod['id']+'] Ellipse angle: '+str(opt_angle)[:6]+'; Ellipse axial ratio: '+str(opt_axial_ratio)[:6]

    # Clean garbage, record results to pod, and return
    gc.collect()
    pod['cutout_clip'] = clip_value
    pod['opt_axial_ratio'] = opt_axial_ratio
    pod['opt_angle'] = opt_angle
    pod['semimaj_initial_pix'] = semimaj_initial
    pod['semimaj_initial_arcsec'] = semimaj_initial * pod['pix_size']
    return pod





# Define function that determines the size of the source aperture in this band
def ApertureSize(pod, verbose=False):
    if pod['verbose']: print '['+pod['id']+'] Commencing determination of appropriate size for source aperture.'
    band_dict = pod['band_dict']
    verbose = pod['verbose']

    # Define sub-function that determines SNR of a defined annulus; and if requested, determines residual between the SNR of a defined annulus, and a target SNR of 2
    def AnnulusSNR(semimaj, pod, cutout, width, i_trans, j_trans, residual):
        semimaj = semimaj[0]
        sig_annulus = ChrisFuncs.Photom.AnnulusQuickSum(cutout, semimaj, width, pod['opt_axial_ratio'], pod['opt_angle'], pod['centre_i'], pod['centre_j'], i_trans, j_trans)
        sig_value = np.median(sig_annulus[2]) #sig_annulus[0] / sig_annulus[1]
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



    # Assuming enough RAM is free, smooth image prior to finding source size (in a way that manages NaNs and replaces them at the end)
    mem_wait = True
    while mem_wait:
        mem_wait = CAAPR_IO.MemCheck(pod['in_fitspath_size'], verbose=verbose)
        if mem_wait:
            print '['+pod['id']+'] Waiting for some RAM to free up before commencing processing.'
            time.sleep(10.0+(5.0*np.random.rand()))
        elif not mem_wait:
            break
    if verbose: print '['+pod['id']+'] Convolving map to lower resolution for radial analysis.'
    pix_size = pod['pix_size']
    res_in = band_dict['beam_width']
    res_out = 3.0*band_dict['beam_width']#36.0
    kernel_fwhm = np.sqrt( (res_out/pix_size)**2.0 - (res_in/pix_size)**2.0 )
    kernel = astropy.convolution.kernels.Gaussian2DKernel(kernel_fwhm)
    #kernel = astropy.convolution.AiryDisk2DKernel(kernel_fwhm)
    cutout_conv = astropy.convolution.convolve_fft(pod['cutout'], kernel, interpolate_nan=True, normalize_kernel=True, ignore_edge_zeros=False)
    cutout_conv[ np.where( np.isnan(pod['cutout'])==True ) ] = np.NaN
    pod['cutout'] = cutout_conv

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
    #[ AnnulusSNR([rad], pod, pod['cutout'], ann_brute_width, i_trans, j_trans, False) for rad in ann_brute_range ]
    snr_success = False
    for i in range(0, len(ann_brute_range)):
        ann_brute_snr = AnnulusSNR([ann_brute_range[i]], pod, pod['cutout'], ann_brute_width, i_trans, j_trans, False)
        if ann_brute_snr>2:
            snr_success = True
            #ann_brute_semimaj = ann_brute_range[i-1]
            ann_bounds = [( ann_brute_range[ max(i-2,0) ], ann_brute_range[ min(i+1,len(ann_brute_range)-1) ] )]
            if verbose: print '['+pod['id']+'] Course analysis finds that radial SNR=2 between semi-major axes of '+str(ann_bounds[0][0]*pod['pix_size'])[:7]+' and '+str(ann_bounds[0][1]*pod['pix_size'])[:7]+' arcseconds.'
            break

    # If SNR=2 threshold not reached, set to minimum semi-major axis of two beam-widths
    if snr_success==False:
        ann_bounds = [( ann_brute_range[i], np.floor(pod['cutout'].shape[0]/2.0) )]
        opt_semimaj_pix = pod['beam_pix'] * 2.0
        opt_semimaj_arcsec = opt_semimaj_pix * pod['pix_size']
        if verbose: print '['+pod['id']+'] No SNR=2 threshold found; hence reverting to two beam-width minimum value of '+str(opt_semimaj_arcsec)[:7]+' arcseconds.'

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
    opt_semimaj_arcsec = opt_semimaj_pix * pod['pix_size']
    if verbose: print '['+pod['id']+'] Precision analysis finds that radial SNR=2 at semi-major axis of '+str(opt_semimaj_arcsec)[:7]+' arcseconds.'

    # For small sources, default to minimum semi-major axis of two beam-widths
    if opt_semimaj_pix<(ann_beams*2.0):
        opt_semimaj_pix = pod['beam_pix'] * 2.0
        opt_semimaj_arcsec = opt_semimaj_pix * pod['pix_size']
        if verbose: print '['+pod['id']+'] Semi-major axis at which SNR=2 is less than two beam-widths; hence reverting to two beam-width minimum value of '+str(opt_semimaj_arcsec)[:7]+' arcseconds.'

    # Deconvolve aperture semi-major axis with beam, by subtracting in quadrature
    adj_semimaj_arcsec = abs( opt_semimaj_arcsec**2.0 - (0.5*band_dict['beam_width'])**2.0 )**0.5
    opt_semimin_arcsec = opt_semimaj_arcsec / pod['opt_axial_ratio']
    adj_semimin_arcsec = abs( opt_semimin_arcsec**2.0 - (0.5*band_dict['beam_width'])**2.0 )**0.5
    adj_ax_ratio = adj_semimaj_arcsec / adj_semimin_arcsec

    # Clean garbage, record results to pod, and return
    gc.collect()
    pod['opt_semimaj_arcsec'] = adj_semimaj_arcsec
    pod['opt_semimaj_pix'] = adj_semimaj_arcsec / pod['pix_size']
    pod['opt_axial_ratio'] = adj_ax_ratio
    return pod





# Define function that combines a set of apertures for a given source into
def CombineAperture(aperture_output_list, source_dict, verbose=False):
    if verbose: print '['+source_dict['name']+'] Combining apertures from all bands to generate final aperture.'



    # Extract various aperture values
    semimaj_arcsec_list = []
    axial_ratio_list = []
    angle_list = []
    for aperture in aperture_output_list:
        semimaj_arcsec_list.append( aperture['opt_semimaj_arcsec'] )
        axial_ratio_list.append( aperture['opt_axial_ratio'] )
        angle_list.append( aperture['opt_angle'] )

    # Find largest semi-major axis, and use to define size of enclosisity array (which will have pixels 1/100th the size of the smallest semi-major axis)
    semimaj_max = np.max(semimaj_arcsec_list)
    semimaj_min = np.min(semimaj_arcsec_list)
    ap_array_pix_size = 0.01 * semimaj_min
    ap_array_scale = int( np.round( semimaj_max / ap_array_pix_size ) )
    ap_array = np.zeros([ 1+(2.2*ap_array_scale), 1+(2.2*ap_array_scale) ])
    centre_i, centre_j = 1+(1.1*ap_array_scale), 1+(1.1*ap_array_scale)
    semimaj_pix_list = np.array(semimaj_arcsec_list) / ap_array_pix_size

    # Loop over each aperture, adding to enclosisity array
    for a in range(0, len(semimaj_pix_list)):
        ap_mask = ChrisFuncs.EllipseMask(ap_array, semimaj_pix_list[a], axial_ratio_list[a], angle_list[a], centre_i, centre_j)
        ap_array[ np.where( ap_mask==1 ) ] += 1
    #ChrisFuncs.Cutout(ap_array, '/home/saruman/spx7cjc/DustPedia/Ap.fits')

    # Fine ellipse that traces edge of enclosisity region
    cont_rad_initial_pix = ( semimaj_min / np.max(axial_ratio_list) ) / ap_array_pix_size
    cont_array = ChrisFuncs.Photom.ContiguousPixels(ap_array, cont_rad_initial_pix, centre_i, centre_j, 0.1)
    cont_x = ((np.where(cont_array==1))[1])
    cont_y = ((np.where(cont_array==1))[0])
    if cont_x.shape[0]>10:
        try:
            cont_ellipse = ChrisFuncs.Photom.EllipseFit(cont_x, cont_y)
            cont_axial_ratio = max(cont_ellipse[1]) / min(cont_ellipse[1])
            cont_angle = cont_ellipse[2]
            cont_semimaj_pix = max([ cont_ellipse[1].max(), semimaj_pix_list.max() ])
        except:
            cont_axial_ratio = 1.0
            cont_angle = 0.0
            cont_semimaj_pix = max([ cont_ellipse[1].max(), semimaj_pix_list.max() ])

    # Convert final semi-major axis back to arcsec, clean garbage, and return results
    cont_semimaj_arcsec = cont_semimaj_pix * ap_array_pix_size
    gc.collect()
    return [cont_semimaj_arcsec, cont_axial_ratio, cont_angle, ap_array]



