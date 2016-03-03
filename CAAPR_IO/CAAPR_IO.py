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
import numpy as np
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



# Dunction that uses band table CSV file to create a dictionary of band values
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



# Function that produces a cutout of a given source in a given band
def Cutout(source_dict, band_dict, output_dir_path, temp_dir_path):

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
        raise ValueError('No appropriately-named input file found in target directroy (please ensure that filesnames are in \"[NAME]_[BAND].fits\" format.')

    # If error maps are being used, construct this path also
    if band_dict['use_error_map']==True:
        in_fitspath_error = in_fitspath+'_Error'
        if os.path.exists(in_fitspath_error+'.fits'):
            in_fitspath_error = in_fitspath_error+'.fits'
        elif os.path.exists(in_fitspath_error+'.fits.gz'):
            in_fitspath_error = in_fitspath_error+'.fits.gz'
        else:
            raise ValueError('No appropriately-named error file found in target directroy (please ensure that error filesnames are in \"[NAME]_[BAND]_Error.fits\" format.')
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
    # Check if Montage is installed
    try:
        import montage_wrapper
        montage_installed = True
    except:
        print 'Montage and/or the Python Montage wrapper module could not be imported. A cruder cutout method will be used, which could cause projection effects for large maps etc. If this matters to you, install Montage and the Python Montage wrapper module!'
        montage_installed = False

    # Construct output path (likewise for error map, if necessary)
    out_fitspath = os.path.join( temp_dir_path, 'Cutouts', source_dict['name'], source_dict['name']+'_'+band_dict['band_name']+'.fits' )
    if band_dict['use_error_map']==True:
        out_fitspath_error = os.path.join( temp_dir_path, 'Cutouts', source_dict['name'], source_dict['name']+'_'+band_dict['band_name']+'_Error.fits' )

    # If Montage is installed, use it to produce cutout using mSubimage function
    if montage_installed:
        montage_wrapper.commands.mSubimage(in_fitspath, out_fitspath, source_dict['ra'], source_dict['dec'], float(band_dict['make_cutout'])/3600.0, hdu=0)
        if band_dict['use_error_map']==True:
            montage_wrapper.commands.mSubimage(in_fitspath_error, out_fitspath_error, source_dict['ra'], source_dict['dec'], float(band_dict['make_cutout'])/3600.0, hdu=0)

    # If montage isn't installed, use the ChrisFuncs cutout function instead
    elif not montage_installed:
        ChrisFuncs.FitsCutout(in_fitspath, source_dict['ra'], source_dict['dec'], int(round(float(band_dict['make_cutout'])/2.0)), exten=0, outfile=out_fitspath)
        if band_dict['use_error_map']==True:
            ChrisFuncs.FitsCutout(in_fitspath_error, source_dict['ra'], source_dict['dec'], int(round(float(band_dict['make_cutout'])/2.0)), exten=0, outfile=out_fitspath_error)

    # Return the directory of the newly-created cutout
    out_fitsdir = os.path.split(out_fitspath)[0]
    return out_fitsdir



# Define function that checks whether a decent amount of RAM is free before allowing things to progress
def MemCheck(pod, thresh_fraction=0.75, thresh_factor=15.0, swap_thresh_fraction=0.2, return_status=False):

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
        if wait_count>=50:
            mem_wait = False

        # If required, conduct wait for a semi-random period of time; otherwise, break
        if mem_wait==True:
            if wait_initial:
                if pod['verbose']: print '['+pod['id']+'] Waiting for necessary RAM to free up before continuing processing.'
                wait_initial = False
            wait_count += 1
            time.sleep(5.0+(10.0*np.random.rand()))
        else:
            break

        # Return memory status
        if return_status:
            return mem_stats



# Define function that makes a grid of operture thumbnails for a given source's output
def ApertureThumbGrid(source_dict, bands_dict, kwargs_dict, aperture_list, aperture_combined):
    import warnings
    warnings.filterwarnings('ignore')

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

    # Begin loop
    for band_name in bands_list:
        for w in range(0, thumb_files):
            if aperture_list[w]['band_name']==band_name:
                b = w
                break
            continue
        thumb_rad = aperture_combined[0] * bands_dict[band_name]['annulus_outer']

        # Assuming cutout exists, appropriately-sized cutoutin current band
        img_input = os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'.fits')
        if not os.path.exists(img_input):
            continue
        img_output = os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'_Thumbnail.fits')
        ChrisFuncs.FitsCutout( img_input, source_dict['ra'], source_dict['dec'], np.ceil(thumb_rad*1.1), exten=0, outfile=img_output )

        # Calculate subplot coordinates
        x_min = x_fig_margin + ( np.mod(float(counter), x_grid) * x_fig_subdim ) + ( np.mod(float(counter), x_grid) * x_fig_margin )
        y_min = 1.0 - ( 1.0 * y_fig_margin ) - ( y_fig_margin + y_fig_subdim + ( np.floor(float(counter)/x_grid) * y_fig_subdim ) + ( np.floor(float(counter)/x_grid) * y_fig_margin ) )
        dx = x_fig_subdim
        dy = y_fig_subdim

        # Disable console output whilst plotting images
        try:
            sys.stdout = open(os.devnull, "w")

            # Create and format image
            vars()['subfig'+str(b)] = aplpy.FITSFigure(img_output, figure=fig, subplot=[x_min, y_min, dx, dy])
            vars()['subfig'+str(b)].show_colorscale(cmap=bands_dict[band_name]['colour_map'], stretch='arcsinh', pmin=7.5, pmax=99.5)
            vars()['subfig'+str(b)].set_nan_color('black')
            vars()['subfig'+str(b)].axis_labels.hide()
            vars()['subfig'+str(b)].tick_labels.hide()
            vars()['subfig'+str(b)].ticks.hide()

            # Plot band-specific aperture
            line_width = 4.0
            band_ap_semimaj = (2.0*aperture_list[b]['opt_semimaj_arcsec'])/3600.0
            band_ap_axial_ratio = aperture_list[b]['opt_axial_ratio']
            band_ap_angle = aperture_list[b]['opt_angle']
            vars()['subfig'+str(b)].show_ellipses(source_dict['ra'], source_dict['dec'], band_ap_semimaj, band_ap_semimaj/band_ap_axial_ratio, angle=band_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width/2.0, linestyle='dotted')

            # Plot combined aperture and sky annulus
            comb_ap_semimaj = (2.0 * aperture_combined[0])/3600.0
            vars()['subfig'+str(b)].show_ellipses(source_dict['ra'], source_dict['dec'], comb_ap_semimaj, comb_ap_semimaj/aperture_combined[1], angle=aperture_combined[2], edgecolor='#00FF40', facecolor='none', linewidth=line_width)
            comb_ann_inner_semimaj = (aperture_combined[0] * 2.0 * bands_dict[band_name]['annulus_inner'])/3600.0
            comb_ann_outer_semimaj = (aperture_combined[0] * 2.0 * bands_dict[band_name]['annulus_outer'])/3600.0
            vars()['subfig'+str(b)].show_ellipses(source_dict['ra'], source_dict['dec'], comb_ann_inner_semimaj, comb_ann_inner_semimaj/aperture_combined[1], angle=aperture_combined[2], edgecolor='#00FF40', facecolor='none', linewidth=line_width/3.0, linestyle='dashed')
            vars()['subfig'+str(b)].show_ellipses(source_dict['ra'], source_dict['dec'], comb_ann_outer_semimaj, comb_ann_outer_semimaj/aperture_combined[1], angle=aperture_combined[2], edgecolor='#00FF40', facecolor='none', linewidth=line_width/3.0, linestyle='dashed')

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
    fig.savefig( os.path.join(kwargs_dict['output_dir_path'],'Aperture_Fitting_Thumbnails',source_dict['name']+'_Thumbnail_Grid.png'), facecolor='white', dpi=75.0)
    pdb.set_trace()
    if kwargs_dict['do_photom']==False:
        shutil.rmtree(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps'))



