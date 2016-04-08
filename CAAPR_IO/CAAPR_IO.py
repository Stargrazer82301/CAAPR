# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import time
import csv
import shutil
import psutil
import multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.io.fits
import astropy.wcs
import aplpy
import ChrisFuncs
plt.ioff()





# Dunction that uses band table CSV file to create a dictionary of band values
def SourcesDictFromCSV(sources_table_path):

    # Initially read in CSV file as a numpy structured array, and prepare output dictionary
    sources_table = np.genfromtxt(sources_table_path, names=True, delimiter=',', dtype=None, comments='#')
    sources_dict = {}

    # Deal with special case of CSV having only 1 line (where looping through lines doesn't work), otherwise do things properly
    if sources_table.shape==():
        sources_dict[ sources_table['name'].tolist() ] = {}
        for field in sources_table.dtype.names:
            sources_dict[ sources_table['name'].tolist() ][field] = sources_table[field].tolist()
    else:

        # Loop over each row of the bands table, using band name as the the key for this 1st-level entry
        for row in range(0, sources_table.shape[0]):
            sources_dict[ sources_table['name'][row] ] = {}

            # Loop over each field in the current row
            for field in sources_table.dtype.names:
                sources_dict[ sources_table['name'][row] ][field] = sources_table[field][row]

    # Return dictionary
    return sources_dict





# Function that uses band table CSV file to create a dictionary of band values
def BandsDictFromCSV(bands_table_path):

    # Initially read in CSV file as a numpy structured array, and prepare output dictionary
    bands_table = np.genfromtxt(bands_table_path, names=True, delimiter=',', dtype=None, comments='#')
    bands_dict = {}

    # Deal with special case of CSV having only 1 line (where looping through lines doesn't work), otherwise do things properly
    if bands_table.shape==():
        bands_dict[ bands_table['band_name'].tolist() ] = {}
        for field in bands_table.dtype.names:
            bands_dict[ bands_table['band_name'].tolist() ][field] = bands_table[field].tolist()
    else:

        # Loop over each row of the bands table, using band name as the the key for this 1st-level entry
        for row in range(0, bands_table.shape[0]):
            bands_dict[ bands_table['band_name'][row] ] = {}

            # Loop over each field in the current row
            for field in bands_table.dtype.names:
                bands_dict[ bands_table['band_name'][row] ][field] = bands_table[field][row]

    # Return dictionary
    return bands_dict





# Function that prepares output file for photometry
def PhotomTablePrepare(bands_dict, kwargs_dict):

    # Use band input table to establish order in which to put bands in header
    bands_table = np.genfromtxt(kwargs_dict['bands_table_path'], delimiter=',', names=True, dtype=None)
    bands_list = bands_table['band_name']

    # Create header, handling special case of a single band
    photom_table_header = 'name'
    if bands_list.shape==():
        bands_list = [bands_list.tolist()]
    for band in bands_list:
        photom_table_header += ','+band+','+band+'_ERR'
    photom_table_header += '\n'

    photom_table_file = open(kwargs_dict['photom_table_path'], 'a')
    photom_table_file.write(photom_table_header)
    photom_table_file.close()





# Function that produces a cutout of a given source in a given band
def Cutout(source_dict, band_dict, output_dir_path, temp_dir_path):
    source_id = source_dict['name']+'_'+band_dict['band_name']

    # Determine whether the user is specificing a directroy full of FITS files in this band (in which case use standardised filename format), or just a single FITS file
    if os.path.isdir(band_dict['band_dir']):
        in_fitspath = os.path.join( band_dict['band_dir'], source_dict['name']+'_'+band_dict['band_name'] )
    elif os.path.isfile(band_dict['band_dir']):
        in_fitspath = os.path.join( band_dict['band_dir'] )

    # Make sure appropriate cutout sub-directories exist in temp directory
    if not os.path.exists( os.path.join( temp_dir_path, 'Cutouts' ) ):
        os.mkdir( os.path.join( temp_dir_path, 'Cutouts' ) )
    if not os.path.exists( os.path.join( temp_dir_path, 'Cutouts', source_dict['name'] ) ):
        os.mkdir( os.path.join( temp_dir_path, 'Cutouts', source_dict['name'] ) )

    # Work out whether the file extension is .fits or .fits.gz
    if os.path.exists(in_fitspath+'.fits'):
        in_fitspath = in_fitspath+'.fits'
    elif os.path.exists(in_fitspath+'.fits.gz'):
        in_fitspath = in_fitspath+'.fits.gz'
    else:
        in_fitspath = None
        print '['+source_id+'] No appropriately-named input file found in target directroy (please ensure that filesnames are in \"[NAME]_[BAND].fits\" format.)'
        print '['+source_id+'] Assuming no data in this band for current source.'
    if in_fitspath!=None:
        print '['+source_id+'] Creating cutout with '+str(band_dict['make_cutout'])+' arcsec diameter.'

        # If error maps are being used, construct this path also
        if band_dict['use_error_map']==True:
            in_fitspath_error = in_fitspath.replace('.fits','_Error.fits')
            if os.path.exists(in_fitspath_error):
                pass
            elif os.path.exists(in_fitspath_error+'.gz'):
                in_fitspath_error = in_fitspath_error+'.gz'
            else:
                raise ValueError('No appropriately-named error file found in target directroy (please ensure that error filesnames are in \"[NAME]_[BAND]_Error.fits\" format.')

        # Determine pixel size
        in_header = astropy.io.fits.getheader(in_fitspath)
        in_wcs = astropy.wcs.WCS(in_header)
        in_cdelt = in_wcs.wcs.cdelt.max()
        in_cdelt_arcsec = in_cdelt * 3600.0

        # Check if Montage is installed
        try:
            import montage_wrapper
            montage_installed = True
        except:
            print '['+source_id+'] Montage and/or the Python Montage wrapper module could not be imported. A cruder cutout method will be used, which could cause projection effects for large maps etc. If this matters to you, install Montage and the Python Montage wrapper module!'
            montage_installed = False

        # Construct output path (likewise for error map, if necessary)
        out_fitspath = os.path.join( temp_dir_path, 'Cutouts', source_dict['name'], source_dict['name']+'_'+band_dict['band_name']+'.fits' )
        if band_dict['use_error_map']==True:
            out_fitspath_error = os.path.join( temp_dir_path, 'Cutouts', source_dict['name'], source_dict['name']+'_'+band_dict['band_name']+'_Error.fits' )

        # If Montage is installed, use it to produce cutout using mSubimage function
        if montage_installed:
            montage_wrapper.commands.mHdr(str(source_dict['ra'])+' '+str(source_dict['dec']), 1.0, os.path.join(temp_dir_path,'hdr.txt'), pix_size=in_cdelt_arcsec)
            #montage_wrapper.commands.mSubimage(in_fitspath, out_fitspath, source_dict['ra'], source_dict['dec'], float(band_dict['make_cutout'])/3600.0, hdu=0)
            montage_wrapper.wrappers.reproject(in_fitspath, out_fitspath, header=os.path.join(temp_dir_path,'hdr.txt'), north_aligned=True, exact_size=True, hdu=0, cleanup=True, silent_cleanup=True)
            if band_dict['use_error_map']==True:
                #montage_wrapper.commands.mSubimage(in_fitspath_error, out_fitspath_error, source_dict['ra'], source_dict['dec'], float(band_dict['make_cutout'])/3600.0, hdu=0)
                montage_wrapper.wrappers.reproject(in_fitspath_error, out_fitspath_error, header=os.path.join(temp_dir_path,'hdr.txt'), north_aligned=True, exact_size=True, hdu=0, cleanup=True, silent_cleanup=True)

        # If montage isn't installed, use the ChrisFuncs cutout function instead
        elif not montage_installed:
            ChrisFuncs.FitsCutout(in_fitspath, source_dict['ra'], source_dict['dec'], int(round(float(band_dict['make_cutout'])/2.0)), exten=0, outfile=out_fitspath)
            if band_dict['use_error_map']==True:
                ChrisFuncs.FitsCutout(in_fitspath_error, source_dict['ra'], source_dict['dec'], int(round(float(band_dict['make_cutout'])/2.0)), exten=0, outfile=out_fitspath_error)

        # Return the directory of the newly-created cutout
        out_fitsdir = os.path.split(out_fitspath)[0]
        return out_fitsdir





# Define function that checks whether a decent amount of RAM is free before allowing things to progress
def MemCheck(pod, thresh_fraction=0.75, thresh_factor=20.0, swap_thresh_fraction=0.5, return_status=False):

    # Start infinite loop
    wait_initial = True
    wait_count = 0
    while True:
        mem_wait = False

        # Assess how much RAM is free
        mem_stats = psutil.virtual_memory()
        mem_usage = 1.0-(float(mem_stats[1])/float(mem_stats[0])) #float(mem_stats[2]) / 100.0

        # Also, require wait if the amount of RAM free is more than some multiple the size of the current file
        mem_free = float(psutil.virtual_memory()[4])
        if thresh_factor!=None:
            if ( float(thresh_factor) * float(pod['in_fitspath_size']) )>mem_free:
                if pod['in_fitspath_size']!=None:
                    if (pod['in_fitspath_size']*thresh_factor)>(0.75):
                        thresh_fraction = 0.25
                    else:
                        mem_wait = True

        # Require wait if less than some fraction of RAM is free
        if thresh_fraction!=None:
            if mem_usage>=float(thresh_fraction):
                mem_wait = True

        # Also, require that some fraction of the swap is free
        swap_stats = psutil.swap_memory()
        swap_usage = float(swap_stats[3]) / 100.0
        if swap_usage>swap_thresh_fraction:
             mem_wait = True

        # If process has waited loads of times, progress regardless
        if wait_count>=20:
            mem_wait = False

        # If required, conduct wait for a semi-random period of time; otherwise, break
        if mem_wait==True:
            if wait_initial:
                if pod['verbose']: print '['+pod['id']+'] Waiting for necessary RAM to free up before continuing processing.'
                wait_initial = False
            wait_count += 1
            time.sleep(10.0+(10.0*np.random.rand()))
        else:
            break

        # Return memory status
        if return_status:
            return mem_stats





# Define function that makes a grid of aperture thumbnails for a given source's output
def ApertureThumbGrid(source_dict, bands_dict, kwargs_dict, aperture_list, aperture_combined):
    import warnings
    warnings.filterwarnings('ignore')
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Producing image grid of aperture thumbnails.'

    # Define sub-function to find all possible factors for given value
    def Factors(value):
        factors = []
        for i in range(1, int(value**0.5)+1):
            if value % i == 0:
                factors.append([i, value/i])
        return factors



    # Find how many thumnails need tiling
    list_files = np.array( os.listdir( os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps') ) )
    thumb_files = sum( [ '.fits' in list_file for list_file in list_files ] )

    # Find grid dimensions that closely match the golden ratio for the number of images being plotted
    phi = (1.0+5.0**0.5)/2.0
    edge_short = np.floor(float(thumb_files)**0.5)
    edge_long = np.round(phi*float(edge_short))
    if int(edge_short*edge_long)<thumb_files:
        edge_long += 1
    if (int(edge_short*edge_long)-int(edge_long))>=thumb_files:
        edge_long -= 1
    if int(edge_long)==3:
        pdb.set_trace()

    # Set up various variables
    counter = 0
    row_counter = 0
    column_counter = 0

    # Calculate figure coordinates (Border of 1, each subplot of width 10, separated by 1)
    x_grid = edge_long
    y_grid = edge_short
    x_subdim = 4.0
    y_subdim = 4.0
    x_margin = 0.4
    y_margin = 0.4
    x_dim = (x_subdim * x_grid) + (x_margin * (x_grid+1.0))
    y_dim = (y_subdim * y_grid) + (y_margin * (y_grid+1.0)) + (1.0 * y_margin)
    x_fig_subdim = x_subdim / x_dim;
    y_fig_subdim = y_subdim / y_dim
    x_fig_margin = x_margin / x_dim
    y_fig_margin = y_margin / y_dim

    # Use band input table to establish order in which to plot thumbnails
    bands_table = np.genfromtxt(kwargs_dict['bands_table_path'], delimiter=',', names=True, dtype=None)
    bands_list = bands_table['band_name']

    # Handle case where there's only one band
    if len(bands_list.shape)==0:
        bands_list = [bands_list.tolist()]

    # Prepare figure and add title
    fig = plt.figure(figsize=(x_dim, y_dim))
    x_title = 1.0 * x_margin #x_fig_margin + ( np.mod(0.0, x_grid) * x_fig_subdim ) + ( np.mod(0.0, x_grid) * x_fig_margin )
    y_title = y_dim - (0.65 * y_margin) #1.0 - ( y_fig_margin + y_fig_subdim + ( np.floor(float(counter)/x_grid) * y_fig_subdim ) + ( np.floor(float(counter)/x_grid) * y_fig_margin ) )
    plt.figtext(x_title/x_dim, y_title/y_dim, source_dict['name'], size=30, color='black', weight='bold', horizontalalignment='left', verticalalignment='top', figure=fig)



    # Begin preliminary loop, to produce cutouts for thumbnails
    thumb_pool = mp.Pool()
    bands_list_present = []
    for band_name in bands_list:

        # Check whether cutout exists,
        img_input = os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'.fits')
        if not os.path.exists(img_input):
            continue
        else:
            bands_list_present.append(band_name)

        # Calculate desired size of thumbnail, based on fitted apertures, and size of map
        thumb_rad = np.ceil( 1.1 * aperture_combined[0] * bands_dict[band_name]['annulus_outer'] )
        img_header = astropy.io.fits.getheader(img_input)
        img_wcs = astropy.wcs.WCS(img_header)
        img_pix_arcsec = np.max(3600.0*img_wcs.wcs.cdelt)
        img_naxis_pix = np.max([img_header['NAXIS1'],img_header['NAXIS1']])
        img_naxis_arcsec = float(img_naxis_pix) * float(img_pix_arcsec)
        img_rad_arcsec = img_naxis_arcsec / 2.0
        if thumb_rad>img_rad_arcsec:
            thumb_rad = img_rad_arcsec

        # Produce cutouts, and end loop
        thumb_pool.apply_async( ThumbCutout, args=(source_dict, bands_dict[band_name], kwargs_dict, img_input, thumb_rad,) )
    thumb_pool.close()
    thumb_pool.join()



    # Begin main thumbnail plotting loop
    aperture_name_list = [ aperture_list[b]['band_name'] for b in range(0, len(aperture_list))  ]
    for band_name in bands_list_present:
        for w in range(0, thumb_files):
            if aperture_list[w]['band_name']==band_name:
                b = w
                break
            continue

        # Check whether cutout exists,
        img_input = os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'.fits')
        if not os.path.exists(img_input):
            continue
        img_output = os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'_Thumbnail.fits')

        # Calculate subplot coordinates
        x_min = x_fig_margin + ( np.mod(float(counter), x_grid) * x_fig_subdim ) + ( np.mod(float(counter), x_grid) * x_fig_margin )
        y_min = 1.0 - ( 1.0 * y_fig_margin ) - ( y_fig_margin + y_fig_subdim + ( np.floor(float(counter)/x_grid) * y_fig_subdim ) + ( np.floor(float(counter)/x_grid) * y_fig_margin ) )
        dx = x_fig_subdim
        dy = y_fig_subdim

        # Disable console output and warnings whilst plotting images
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sys.stdout = open(os.devnull, "w")

                # Create and format image
                vars()['subfig'+str(b)] = aplpy.FITSFigure(img_output, figure=fig, subplot=[x_min, y_min, dx, dy])
                vars()['subfig'+str(b)].show_colorscale(cmap=bands_dict[band_name]['colour_map'], stretch='arcsinh', pmin=7.5, pmax=99.5)
                vars()['subfig'+str(b)].set_nan_color('black')
                vars()['subfig'+str(b)].axis_labels.hide()
                vars()['subfig'+str(b)].tick_labels.hide()
                vars()['subfig'+str(b)].ticks.hide()

                # Extract band-specific aperture dimensions
                band_ap_angle = aperture_list[b]['opt_angle']
                band_ap_axial_ratio = aperture_list[b]['opt_axial_ratio']
                band_ap_semimaj = (aperture_list[b]['opt_semimaj_arcsec'])/3600.0
                band_ap_semimin = band_ap_semimaj / band_ap_axial_ratio
                band_beam_width = bands_dict[band_name]['beam_arcsec'] / 3600.0
                print band_name
                print band_ap_semimaj*3600.0
                band_ap_semimaj = ( band_ap_semimaj**2.0 + (0.5*band_beam_width)**2.0 )**0.5
                print band_ap_semimaj*3600.0
                band_ap_semimin = ( band_ap_semimin**2.0 + (0.5*band_beam_width)**2.0 )**0.5
                band_ap_axial_ratio = band_ap_semimaj / band_ap_semimaj

                # Plot band-specific aperture (if one was provided)
                line_width = 4.0
                if bands_dict[band_name]['consider_aperture']==True:
                    vars()['subfig'+str(b)].show_ellipses(source_dict['ra'], source_dict['dec'], 2.0*band_ap_semimaj, 2.0*band_ap_semimin, angle=band_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width/2.0, linestyle='dotted')

                # Extract combined aperture dimensions
                comb_ap_angle = aperture_combined[2]
                comb_ap_axial_ratio = aperture_combined[1]
                comb_ap_semimaj = aperture_combined[0]/3600.0
                comb_ap_semimin = comb_ap_semimaj / comb_ap_axial_ratio
                print comb_ap_semimaj*3600.0
                comb_ap_semimaj = ( comb_ap_semimaj**2.0 + band_beam_width**2.0 )**0.5
                print comb_ap_semimaj*3600.0
                comb_ap_semimin = ( comb_ap_semimin**2.0 + band_beam_width**2.0 )**0.5

                # Plot combined aperture
                vars()['subfig'+str(b)].show_ellipses(source_dict['ra'], source_dict['dec'], 2.0*comb_ap_semimaj, 2.0*comb_ap_semimin, angle=comb_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width)

                # Plot combined background annulus
                band_ann_inner_semimaj = comb_ap_semimaj * 2.0 * bands_dict[band_name]['annulus_inner']
                band_ann_outer_semimaj = comb_ap_semimaj * 2.0 * bands_dict[band_name]['annulus_outer']
                vars()['subfig'+str(b)].show_ellipses(source_dict['ra'], source_dict['dec'], band_ann_inner_semimaj, band_ann_inner_semimaj/comb_ap_axial_ratio, angle=comb_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width/3.0, linestyle='dashed')
                vars()['subfig'+str(b)].show_ellipses(source_dict['ra'], source_dict['dec'], band_ann_outer_semimaj, band_ann_outer_semimaj/comb_ap_axial_ratio, angle=comb_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width/3.0, linestyle='dashed')

                # Plot label
                vars()['subfig'+str(b)].add_label(0.035, 0.92, bands_dict[band_name]['band_name'], relative=True, size=20, color='white', horizontalalignment='left')

                # Restore console output
                sys.stdout = sys.__stdout__
        except:
            sys.stdout = sys.__stdout__
            pdb.set_trace()

        # Progress counters
        counter += 1
        column_counter += 1
        if np.mod(float(counter)+1, x_grid)==0:
            row_counter += 1
            column_counter = 0

    # Save figure, and remove temporary files
    fig.savefig( os.path.join(kwargs_dict['output_dir_path'],'Aperture_Fitting_Thumbnails',source_dict['name']+'_Thumbnail_Grid.png'), facecolor='white', dpi=100.0)
    [os.remove(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',processed_map)) for processed_map in os.listdir(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps')) if '.fits' in processed_map]
    gc.collect()





# Define function that makes a grid of photometry thumbnails for a given source's output
def PhotomThumbGrid(source_dict, bands_dict, kwargs_dict):
    import warnings
    warnings.filterwarnings('ignore')
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Producing image grid of photometry thumbnails.'

    # Define sub-function to find all possible factors for given value
    def Factors(value):
        factors = []
        for i in range(1, int(value**0.5)+1):
            if value % i == 0:
                factors.append([i, value/i])
        return factors



    # Read in aperture file
    aperture_table = np.genfromtxt(kwargs_dict['aperture_table_path'], delimiter=',', names=True, dtype=None)
    aperture_index = np.where( aperture_table['name']==source_dict['name'] )
    if aperture_index[0].shape[0]>1:
        raise ValueError('Aperture value caontains more than one entry for current galaxy')
    else:
        aperture_index = aperture_index[0][0]

    # Extract aperture corresponding to current source, dealing with special case where aperture file contains only one source
    if np.where( aperture_table['name']==source_dict['name'] )[0][0]==0:
        opt_semimaj_arcsec = aperture_table['semimaj_arcsec'].tolist()
        opt_axial_ratio = aperture_table['axial_ratio'].tolist()
        opt_angle = aperture_table['pos_angle'].tolist()
    else:
        opt_semimaj_arcsec = aperture_table['semimaj_arcsec'][aperture_index]
        opt_axial_ratio = aperture_table['axial_ratio'][aperture_index]
        opt_angle = aperture_table['pos_angle'][aperture_index]



    # Find how many thumnails need tiling
    list_files = np.array( os.listdir( os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps') ) )
    thumb_files = sum( [ '.fits' in list_file for list_file in list_files ] )

    # Find grid dimensions that closely match the golden ratio for the number of images being plotted
    phi = (1.0+5.0**0.5)/2.0
    edge_short = np.floor(float(thumb_files)**0.5)
    edge_long = np.round(phi*float(edge_short))
    if int(edge_short*edge_long)<thumb_files:
        edge_long += 1
    if (int(edge_short*edge_long)-int(edge_long))>=thumb_files:
        edge_long -= 1

    # Set up various variables
    counter = 0
    row_counter = 0
    column_counter = 0

    # Calculate figure coordinates (Border of 1, each subplot of width 10, separated by 1)
    x_grid = edge_long
    y_grid = edge_short
    x_subdim = 4.0
    y_subdim = 4.0
    x_margin = 0.4
    y_margin = 0.4
    x_dim = (x_subdim * x_grid) + (x_margin * (x_grid+1.0))
    y_dim = (y_subdim * y_grid) + (y_margin * (y_grid+1.0)) + (1.0 * y_margin)
    x_fig_subdim = x_subdim / x_dim;
    y_fig_subdim = y_subdim / y_dim
    x_fig_margin = x_margin / x_dim
    y_fig_margin = y_margin / y_dim

    # Use band input table to establish order in which to plot thumbnails
    bands_table = np.genfromtxt(kwargs_dict['bands_table_path'], delimiter=',', names=True, dtype=None)
    bands_list = bands_table['band_name']

    # Handle case where there's only one band
    if len(bands_list.shape)==0:
        bands_list = [bands_list.tolist()]

    # Prepare figure and add title
    fig = plt.figure(figsize=(x_dim, y_dim))
    x_title = 1.0 * x_margin #x_fig_margin + ( np.mod(0.0, x_grid) * x_fig_subdim ) + ( np.mod(0.0, x_grid) * x_fig_margin )
    y_title = y_dim - (0.65 * y_margin) #1.0 - ( y_fig_margin + y_fig_subdim + ( np.floor(float(counter)/x_grid) * y_fig_subdim ) + ( np.floor(float(counter)/x_grid) * y_fig_margin ) )
    plt.figtext(x_title/x_dim, y_title/y_dim, source_dict['name'], size=30, color='black', weight='bold', horizontalalignment='left', verticalalignment='top', figure=fig)



    # Begin preliminary loop, to produce cutouts for thumbnails
    thumb_pool = mp.Pool()
    bands_list_present = []
    for band_name in bands_list:

        # Check whether cutout exists,
        img_input = os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'.fits')
        if not os.path.exists(img_input):
            continue
        else:
            bands_list_present.append(band_name)

        # Calculate desired size of thumbnail, based on fitted apertures, and size of map
        thumb_rad = np.ceil( 1.1 * opt_semimaj_arcsec * bands_dict[band_name]['annulus_outer'] )
        img_header = astropy.io.fits.getheader(img_input)
        img_wcs = astropy.wcs.WCS(img_header)
        img_pix_arcsec = np.max(3600.0*img_wcs.wcs.cdelt)
        img_naxis_pix = np.max([img_header['NAXIS1'],img_header['NAXIS1']])
        img_naxis_arcsec = float(img_naxis_pix) * float(img_pix_arcsec)
        img_rad_arcsec = img_naxis_arcsec / 2.0
        if thumb_rad>img_rad_arcsec:
            thumb_rad = img_rad_arcsec

        # Produce cutouts, and end loop
        thumb_pool.apply_async( ThumbCutout, args=(source_dict, bands_dict[band_name], kwargs_dict, img_input, thumb_rad,) )
    thumb_pool.close()
    thumb_pool.join()



    # Begin main thumbnail plotting loop
    w = 0
    for band_name in bands_list_present:
        w += 1

        # Check whether cutout exists,
        img_input = os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'.fits')
        if not os.path.exists(img_input):
            continue
        img_output = os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'_Thumbnail.fits')

        # Calculate subplot coordinates
        x_min = x_fig_margin + ( np.mod(float(counter), x_grid) * x_fig_subdim ) + ( np.mod(float(counter), x_grid) * x_fig_margin )
        y_min = 1.0 - ( 1.0 * y_fig_margin ) - ( y_fig_margin + y_fig_subdim + ( np.floor(float(counter)/x_grid) * y_fig_subdim ) + ( np.floor(float(counter)/x_grid) * y_fig_margin ) )
        dx = x_fig_subdim
        dy = y_fig_subdim

        # Disable console output and warnings whilst plotting images
#        try:
#            with warnings.catch_warnings():
#                warnings.simplefilter("ignore")
#                sys.stdout = open(os.devnull, "w")

        # Create and format image
        vars()['subfig'+str(w)] = aplpy.FITSFigure(img_output, figure=fig, subplot=[x_min, y_min, dx, dy])
        vars()['subfig'+str(w)].show_colorscale(cmap=bands_dict[band_name]['colour_map'], stretch='arcsinh', pmin=7.5, pmax=99.5)
        vars()['subfig'+str(w)].set_nan_color('black')
        vars()['subfig'+str(w)].axis_labels.hide()
        vars()['subfig'+str(w)].tick_labels.hide()
        vars()['subfig'+str(w)].ticks.hide()

        # Extract combined aperture dimensions
        band_beam_width = bands_dict[band_name]['beam_arcsec'] / 3600.0
        comb_ap_angle = opt_angle
        comb_ap_axial_ratio = opt_axial_ratio
        comb_ap_semimaj = opt_semimaj_arcsec/3600.0
        comb_ap_semimin = comb_ap_semimaj / comb_ap_axial_ratio
        comb_ap_semimaj = ( comb_ap_semimaj**2.0 + band_beam_width**2.0 )**0.5
        comb_ap_semimin = ( comb_ap_semimin**2.0 + band_beam_width**2.0 )**0.5

        # Plot combined aperture
        line_width = 4.0
        vars()['subfig'+str(w)].show_ellipses(source_dict['ra'], source_dict['dec'], 2.0*comb_ap_semimaj, 2.0*comb_ap_semimin, angle=comb_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width)

        # Plot combined background annulus
        band_ann_inner_semimaj = comb_ap_semimaj * 2.0 * bands_dict[band_name]['annulus_inner']
        band_ann_outer_semimaj = comb_ap_semimaj * 2.0 * bands_dict[band_name]['annulus_outer']
        vars()['subfig'+str(w)].show_ellipses(source_dict['ra'], source_dict['dec'], band_ann_inner_semimaj, band_ann_inner_semimaj/comb_ap_axial_ratio, angle=comb_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width/3.0, linestyle='dashed')
        vars()['subfig'+str(w)].show_ellipses(source_dict['ra'], source_dict['dec'], band_ann_outer_semimaj, band_ann_outer_semimaj/comb_ap_axial_ratio, angle=comb_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width/3.0, linestyle='dashed')

        # Plot label
        vars()['subfig'+str(w)].add_label(0.035, 0.92, bands_dict[band_name]['band_name'], relative=True, size=20, color='white', horizontalalignment='left')

#                # Restore console output
#                sys.stdout = sys.__stdout__
#        except:
#            sys.stdout = sys.__stdout__
#            pdb.set_trace()

        # Progress counters
        counter += 1
        column_counter += 1
        if np.mod(float(counter)+1, x_grid)==0:
            row_counter += 1
            column_counter = 0

    # Save figure, and remove temporary files
    fig.savefig( os.path.join(kwargs_dict['output_dir_path'],'Photometry_Thumbnails',source_dict['name']+'_Thumbnail_Grid.png'), facecolor='white', dpi=100.0)
    [os.remove(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',processed_map)) for processed_map in os.listdir(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps')) if '.fits' in processed_map]
    gc.collect()





# Define function for producing an appropriately-sized cutout in current band (with console output disabled), with an equinox header keyword to make APLpy *shut up*
def ThumbCutout(source_dict, band_dict, kwargs_dict, img_input, thumb_rad):
    img_output = os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_dict['band_name']+'_Thumbnail.fits')
    try:
        sys.stdout = open(os.devnull, "w")
        ChrisFuncs.FitsCutout( img_input, source_dict['ra'], source_dict['dec'], thumb_rad, exten=0, outfile=img_output )
        cutout_data, cutout_header = astropy.io.fits.getdata(img_output, header=True)
        del cutout_header['RADESYS']
        cutout_header['EQUINOX'] = 2000.0
        cutout_header['EPOCH'] = 2000.00
        astropy.io.fits.writeto(img_output, cutout_data, header=cutout_header, clobber=True)
        sys.stdout = sys.__stdout__
    except:
        sys.stdout = sys.__stdout__
        pdb.set_trace()





"""
# Open FITS file in question (it is assumed that only the zeroth extension is of any interest)
in_fitsdata = astropy.io.fits.open(in_fitspath)
in_image = in_fitsdata[0].data
in_header = in_fitsdata[0].header
in_fitsdata.close()
in_wcs = astropy.wcs.WCS(in_header)
pix_size = in_wcs.wcs.cdelt

# If error maps are being used, open this file also
if band_dict['use_error_map']==True:
    in_fitsdata_error = astropy.io.fits.open(in_fitspath_error)
    in_image_error = in_fitsdata_error[0].data
    in_header_error = in_fitsdata_error[0].header
    in_fitsdata_error.close()
    in_wcs_error = astropy.wcs.WCS(in_header_error)
    pix_size_error = in_wcs_error.wcs.cdelt
"""