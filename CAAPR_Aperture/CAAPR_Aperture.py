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
def PipelineAperture(source_dict, band_dict, kwargs_dict):
    source_id = source_dict['name']+'_'+band_dict['band_name']



    # Determine whether the user is specificing a directroy full of FITS files in this band (in which case use standardised filename format), or just a single FITS file
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
        #raise ValueError('No appropriately-named input file found in target directroy (please ensure that filesnames are in \"[NAME]_[BAND].fits\" format.')
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
               'output_dir_path':kwargs_dict['output_dir_path'],
               'temp_dir_path':kwargs_dict['temp_dir_path'],
               'in_fitspath_size':in_fitspath_size,
               'id':source_id,
               'verbose':kwargs_dict['verbose']}



        # Run pod through preliminary processing, to determine initial quantities; if target not within bounds of map, end processing here
        CAAPR_IO.MemCheck(pod)
        pod = CAAPR_Pipeline.MapPrelim(pod)
        if pod['within_bounds']==False:
            return pod



        # If star-removal is required, run pod through AstroMagic
        if band_dict['remove_stars']==True:
            CAAPR_IO.MemCheck(pod)
            CAAPR_Pipeline.AstroMagic(pod)



        # Run pod through function that determines aperture shape, to provide preliminary estimate to facilitate removal of large-scale sky
        CAAPR_IO.MemCheck(pod)
        pod = ApertureShape(pod)



        # Run pod through function that removes large-scale sky using a 2-dimensional polynomial filter
        CAAPR_IO.MemCheck(pod)
        pod = CAAPR_Pipeline.PolySub(pod)



        # If sky polynomial removed, run pod through function that determines aperture shape, to provide final estiamte
        if pod['sky_poly']!=False:
            CAAPR_IO.MemCheck(pod)
            pod = ApertureShape(pod)



        # Run pod through function that determines aperture size
        CAAPR_IO.MemCheck(pod)
        pod = ApertureSize(pod)



        # If thumbnail images have been requested, save a copy of the current image (ie, with any star and/or background subtaction)
        if kwargs_dict['thumbnails']==True or kwargs_dict['do_photom']==True:
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
    if verbose: print '['+pod['id']+'] Ellipse angle: '+str(opt_angle)[:6]+'; Ellipse axial ratio: '+str(opt_axial_ratio)[:6]+'.'

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



    # Construct kernel with FWHM equal to 3 beam-widths, by which to smooth map
    if verbose: print '['+pod['id']+'] Convolving map to lower resolution for radial analysis.'
    pix_size = pod['pix_size']
    res_in = band_dict['beam_width']
    res_out = 2.0*band_dict['beam_width']#36.0
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
        CAAPR_IO.MemCheck(pod, thresh_fraction=0.66, thresh_factor=20.0)
        kernel = astropy.convolution.kernels.Gaussian2DKernel(kernel_fwhm)
        #kernel = astropy.convolution.AiryDisk2DKernel(kernel_fwhm)
        pod['cutout'] = astropy.convolution.convolve_fft(pod['cutout'], kernel, interpolate_nan=False, normalize_kernel=True, ignore_edge_zeros=False)
        pod['cutout'][ np.where( np.isnan(cutout_unconv)==True ) ] = np.NaN



    # Prepare arrays of transposed coordinates, to allow for rapid radial evaluating
    CAAPR_IO.MemCheck(pod)
    if verbose: print '['+pod['id']+'] Constructing arrays of transposed radial coordinates.'
    coords_trans = ChrisFuncs.Photom.AnnulusQuickPrepare(pod['cutout'], pod['opt_angle'], pod['centre_i'], pod['centre_j'])
    i_trans, j_trans = coords_trans[0], coords_trans[1]

    # To start with, to make new estimate of map noise that isn't contaminated by the target galaxy, by masking all pixels beyond semi-major axis suggested by contiguous significant pixels.
    CAAPR_IO.MemCheck(pod)
    brute_mask = ChrisFuncs.Photom.EllipseMask(pod['cutout'], pod['semimaj_initial_pix'], pod['opt_axial_ratio'], pod['opt_angle'], pod['centre_i'], pod['centre_j'])
    cutout_brute_masked = pod['cutout'].copy()
    cutout_brute_masked[ np.where( brute_mask==1 ) ] = np.nan
    cutout_clip_masked = ChrisFuncs.SigmaClip(cutout_brute_masked, tolerance=0.001, sigma_thresh=3.0, median=True)
    pod['cutout_clip'] = cutout_clip_masked

    # Now, perform a coarse brute force ckeck of a small number of radii over a wide range, to find rough location of edge of the source
    CAAPR_IO.MemCheck(pod)
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
            if verbose: print '['+pod['id']+'] Course analysis finds that radial SNR=2 between semi-major axes of '+str(ann_bounds[0][0]*pod['pix_size'])[:7]+' and '+str(ann_bounds[0][1]*pod['pix_size'])[:7]+' arcseconds.'
            break

    # If SNR=2 threshold not reached, set to minimum semi-major axis of two beam-widths
    if snr_success==False:
        ann_bounds = [( ann_brute_range[i], np.floor(pod['cutout'].shape[0]/2.0) )]
        opt_semimaj_pix = pod['beam_pix'] * 2.0
        opt_semimaj_arcsec = opt_semimaj_pix * pod['pix_size']
        if verbose: print '['+pod['id']+'] No SNR=2 threshold found; hence reverting to two beam-width minimum value of '+str(opt_semimaj_arcsec)[:7]+' arcseconds.'

    # Now use scipy differential evolution optimisation to find, with more precision, the semi-major axis at which annulus falls to a SNR of 2
    CAAPR_IO.MemCheck(pod)
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

    # Establish what fraction of the pixels inside a band's aperture are NaNs
    pix_good = ChrisFuncs.Photom.EllipseQuickSum(pod['cutout'], opt_semimaj_pix, pod['opt_axial_ratio'], pod['opt_angle'], pod['centre_i'], pod['centre_j'], i_trans, j_trans)[1]
    pix_tot = np.where( ChrisFuncs.Photom.EllipseMask(pod['cutout'], opt_semimaj_pix, pod['opt_axial_ratio'], pod['opt_angle'], pod['centre_i'], pod['centre_j']) == 1 )[0].shape[0]
    pix_good_frac = float(pix_good) / float(pix_tot)

    # Before final reporting, tidy up and return to unconvolved cutout
    pod['cutout'] = cutout_unconv
    del(cutout_unconv)
    gc.collect()

    # If more than 10% of the pixels in the aperture are NaNs, report NaN values for aperture dimensions; else proceed normally
    if pix_good_frac<0.9:
        if verbose: print '['+pod['id']+'] More than 10% of pixels in fitted aperture are NaNs; aperture dimensions in this band will not be considered for final aperture.'
        pod['opt_semimaj_arcsec'] = np.NaN
        pod['opt_semimaj_pix'] = np.NaN
        pod['opt_axial_ratio'] = np.NaN
    else:

        # Deconvolve aperture semi-major axis with beam, by subtracting in quadrature
        adj_semimaj_arcsec = abs( opt_semimaj_arcsec**2.0 - (0.5*band_dict['beam_width'])**2.0 )**0.5
        opt_semimin_arcsec = opt_semimaj_arcsec / pod['opt_axial_ratio']
        adj_semimin_arcsec = abs( opt_semimin_arcsec**2.0 - (0.5*band_dict['beam_width'])**2.0 )**0.5
        adj_ax_ratio = adj_semimaj_arcsec / adj_semimin_arcsec

        # Record final dimensions to pod, and return
        pod['opt_semimaj_arcsec'] = adj_semimaj_arcsec
        pod['opt_semimaj_pix'] = adj_semimaj_arcsec / pod['pix_size']
        pod['opt_axial_ratio'] = adj_ax_ratio
        return pod





# Define function that combines a set of apertures for a given source into
def CombineAperture(aperture_output_list, source_dict, kwargs_dict):
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Combining individual apertures from all bands to generate final aperture.'



    # Extract various aperture values
    semimaj_arcsec_list = []
    axial_ratio_list = []
    angle_list = []
    for aperture in aperture_output_list:
        semimaj_arcsec_list.append( aperture['opt_semimaj_arcsec'] )
        axial_ratio_list.append( aperture['opt_axial_ratio'] )
        angle_list.append( aperture['opt_angle'] )

    # Find largest semi-major axis, and use to define size of enclosisity array (which will have pixels some fraction the size of the smallest semi-major axis)
    semimaj_max = np.nanmax(semimaj_arcsec_list)
    semimaj_min = np.nanmin(semimaj_arcsec_list)
    ap_array_pix_size = 0.025 * semimaj_min
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

    # Fine ellipse that traces edge of enclosisity region
    cont_rad_initial_pix = ( semimaj_min / np.nanmax(axial_ratio_list) ) / ap_array_pix_size
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

    # Convert final semi-major axis back to arcsec and apply expanson factor, then clean garbage and return results
    if isinstance(kwargs_dict['expansion_factor'], float) or isinstance(kwargs_dict['expansion_factor'], int):
        expansion_facor = kwargs_dict['expansion_factor']
    else:
        expansion_facor = 1.0
    cont_semimaj_arcsec = cont_semimaj_pix * ap_array_pix_size * expansion_facor

    # Clean garbage and return results
    gc.collect()
    return [cont_semimaj_arcsec, cont_axial_ratio, cont_angle, ap_array]



