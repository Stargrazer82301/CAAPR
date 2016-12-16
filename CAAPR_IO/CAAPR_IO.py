# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import time
import csv
import shutil
import multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy
astropy.log.setLevel('ERROR')
import astropy.io.fits
import astropy.wcs
import aplpy
import ChrisFuncs
import CAAPR
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
                if field=='name':
                    sources_dict[ sources_table['name'][row] ][field] = str(sources_table[field][row])

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





# Function that prepares output directory
def OutputDirPrepare(kwargs_dict):



    # If needed, create output directory (warning user if there is an existing directory there, which hence may have its contents overwritten)
    if os.path.exists( kwargs_dict['output_dir_path'] ):
        print '[CAAPR] Warning: Output directory already exists; some files may be overridden'
    else:
        os.mkdir( kwargs_dict['output_dir_path'] )

    # If aperture fitting requested, make corresponding output sub-directory
    if kwargs_dict['fit_apertures'] and kwargs_dict['thumbnails'] and not os.path.exists( os.path.join( kwargs_dict['output_dir_path'], 'Aperture_Fitting_Thumbnails' ) ):
        os.mkdir( os.path.join(kwargs_dict['output_dir_path'],'Aperture_Fitting_Thumbnails') )

    # If actual photometry requested, make corresponding output sub-directory
    if kwargs_dict['do_photom'] and kwargs_dict['thumbnails'] and not os.path.exists( os.path.join( kwargs_dict['output_dir_path'], 'Photometry_Thumbnails' ) ):
        os.mkdir( os.path.join( kwargs_dict['output_dir_path'],'Photometry_Thumbnails' ) )




# Function that prepares temporary directory
def TempDirPrepare(kwargs_dict):



    # Create temporary directory; if temporary directory already exists, delete it and make a new one
    if os.path.exists( kwargs_dict['temp_dir_path'] ):
        shutil.rmtree( kwargs_dict['temp_dir_path'] )
    os.mkdir( kwargs_dict['temp_dir_path'] )

    # If star-subtraction requested, make temporary sub-directory to hold AstroMagic products
    if kwargs_dict['starsub']==True:
        os.mkdir( os.path.join( kwargs_dict['temp_dir_path'], 'AstroMagic' ) )

    # If thumbnails requested, make temporary sub-directory for thumbnail cutouts
    if kwargs_dict['thumbnails']==True:
        os.mkdir( os.path.join( kwargs_dict['temp_dir_path'],'Processed_Maps' ) )





# Function that prepares output file for apertures
def ApertureTablePrepare(kwargs_dict):



    # Check whether a new aperture file is required; if not, immediately return kwargs dict unchanged
    if kwargs_dict['fit_apertures']!=True:
        return kwargs_dict

    # If no aperture table path is specified construct defualt; otherwis, use path supplied by user
    if kwargs_dict['aperture_table_path']==None:
        kwargs_dict['aperture_table_path'] = os.path.join( kwargs_dict['output_dir_path'], 'CAAPR_Aperture_Table_'+kwargs_dict['timestamp']+'.csv' )

        # Create output file, and write standard header
        aperture_table_header = 'name,semimaj_arcsec,axial_ratio,pos_angle\n'
        aperture_table_file = open( kwargs_dict['aperture_table_path'], 'a')
        aperture_table_file.write(aperture_table_header)
        aperture_table_file.close()

    # Return (potentially-updated) kwargs dict
    return kwargs_dict





# Function that prepares output file for photometry
def PhotomTablePrepare(kwargs_dict):



    # Check that photometry output table is required; if not, immediately return kwargs dict unchanged
    if kwargs_dict['do_photom']!=True:
        return kwargs_dict

    # If no photometry table path is specified construct defualt; otherwis, use path supplied by user
    if kwargs_dict['photom_table_path']==None:
        kwargs_dict['photom_table_path'] = os.path.join( kwargs_dict['output_dir_path'], 'CAAPR_Photom_Table_'+kwargs_dict['timestamp']+'.csv' )

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

    # Create output file, and write output header to it
    photom_table_file = open(kwargs_dict['photom_table_path'], 'a')
    photom_table_file.write(photom_table_header)
    photom_table_file.close()

    # Return (potentially-updated) kwargs dict
    return kwargs_dict





# Function that produces a cutout of a given source in a given band
def Cutout(source_dict, band_dict, kwargs_dict):



    # Check whether cutout has been requested; if not, return band path unchanged
    if str(band_dict['make_cutout'])==True:
        raise Exception('If you want to produce a cutout, please set the \'make_cutout\' field of the band table to be your desired cutout width, in arcsec.')
    if not float(band_dict['make_cutout'])>0:
        return band_dict



    # Determine whether the user is specificing a directroy full of FITS files in this band (in which case use standardised filename format), or just a single FITS file
    if os.path.isdir(band_dict['band_dir']):
        in_fitspath = os.path.join( band_dict['band_dir'], source_dict['name']+'_'+band_dict['band_name'] )
    elif os.path.isfile(band_dict['band_dir']):
        in_fitspath = os.path.join( band_dict['band_dir'] )

    # Make sure appropriate cutout sub-directories exist in temp directory
    if not os.path.exists( os.path.join( kwargs_dict['temp_dir_path'], 'Cutouts' ) ):
        os.mkdir( os.path.join( kwargs_dict['temp_dir_path'], 'Cutouts' ) )
    if not os.path.exists( os.path.join( kwargs_dict['temp_dir_path'], 'Cutouts', source_dict['name'] ) ):
        os.mkdir( os.path.join( kwargs_dict['temp_dir_path'], 'Cutouts', source_dict['name'] ) )

    # Work out whether the file extension is .fits or .fits.gz
    if os.path.exists(in_fitspath+'.fits'):
        in_fitspath = in_fitspath+'.fits'
    elif os.path.exists(in_fitspath+'.fits.gz'):
        in_fitspath = in_fitspath+'.fits.gz'
    else:
        in_fitspath = None
        return band_dict

    # If error maps are being used, construct this path also
    if band_dict['use_error_map']==True:
        in_fitspath_error = in_fitspath.replace('.fits','_Error.fits')
        if os.path.exists(in_fitspath_error):
            pass
        elif os.path.exists(in_fitspath_error+'.gz'):
            in_fitspath_error = in_fitspath_error+'.gz'
        else:
            raise Exception('No appropriately-named error file found in target directroy (please ensure that error filesnames are in \"[NAME]_[BAND]_Error.fits\" format.')

    # Construct output path (likewise for error map, if necessary)
    out_fitspath = os.path.join( kwargs_dict['temp_dir_path'], 'Cutouts', source_dict['name'], source_dict['name']+'_'+band_dict['band_name']+'.fits' )
    if band_dict['use_error_map']==True:
        out_fitspath_error = os.path.join( kwargs_dict['temp_dir_path'], 'Cutouts', source_dict['name'], source_dict['name']+'_'+band_dict['band_name']+'_Error.fits' )

    # Create cutout
    ChrisFuncs.FitsCutout(in_fitspath, source_dict['ra'], source_dict['dec'], int(round(float(band_dict['make_cutout'])/2.0)), exten=0, outfile=out_fitspath, reproj=True)
    if band_dict['use_error_map']==True:
        ChrisFuncs.FitsCutout(in_fitspath_error, source_dict['ra'], source_dict['dec'], int(round(float(band_dict['make_cutout'])/2.0)), exten=0, outfile=out_fitspath_error, reproj=True)

    # Return the directory of the newly-created cutout
    out_fitsdir = os.path.split(out_fitspath)[0]
    band_dict['band_dir'] = out_fitsdir
    return band_dict





# Function that determines if there is unncessary 'padding' around the coverage region of a map, that can be removed
def UnpaddingCutout(source_dict, band_dict, kwargs_dict):



    # Only proceed with unpadding if the user hasn't requested a particular cutout; if they have, return band dict unchanged
    if band_dict['make_cutout']!=False:
        return band_dict



    # Make sure that band directory isn't stuck on cutout directory from previous source
    if 'band_dir_inviolate' in band_dict.keys():
        band_dict['band_dir'] = band_dict['band_dir_inviolate']

    # Determine whether the user is specificing a directroy full of FITS files in this band (in which case use standardised filename format), or just a single FITS file
    if os.path.isdir(band_dict['band_dir']):
        in_fitspath = os.path.join( band_dict['band_dir'], source_dict['name']+'_'+band_dict['band_name'] )
    elif os.path.isfile(band_dict['band_dir']):
        in_fitspath = os.path.join( band_dict['band_dir'] )
    else:
        pdb.set_trace()

    # Make sure appropriate cutout sub-directories exist in temp directory
    if not os.path.exists( os.path.join( kwargs_dict['temp_dir_path'], 'Cutouts' ) ):
        os.mkdir( os.path.join( kwargs_dict['temp_dir_path'], 'Cutouts' ) )
    if not os.path.exists( os.path.join( kwargs_dict['temp_dir_path'], 'Cutouts', source_dict['name'] ) ):
        os.mkdir( os.path.join( kwargs_dict['temp_dir_path'], 'Cutouts', source_dict['name'] ) )

    # Work out whether the file extension is .fits or .fits.gz
    if os.path.exists(in_fitspath+'.fits'):
        in_fitspath = in_fitspath+'.fits'
    elif os.path.exists(in_fitspath+'.fits.gz'):
        in_fitspath = in_fitspath+'.fits.gz'
    else:
        in_fitspath = None
        return band_dict

    # If error maps are being used, construct this path also
    if band_dict['use_error_map']==True:
        in_fitspath_error = in_fitspath.replace('.fits','_Error.fits')
        if os.path.exists(in_fitspath_error):
            pass
        elif os.path.exists(in_fitspath_error+'.gz'):
            in_fitspath_error = in_fitspath_error+'.gz'
        else:
            raise Exception('No appropriately-named error file found in target directroy (please ensure that error filesnames are in \"[NAME]_[BAND]_Error.fits\" format.')

    # Load in map and extract WCS
    in_fits, in_header = astropy.io.fits.getdata(in_fitspath, header=True)
    in_wcs = astropy.wcs.WCS(in_header)

    # Measure size of padding borders
    in_borders = np.where( np.isnan( in_fits )==False )
    x_min_border = np.min(in_borders[1])
    x_max_border = in_fits.shape[1] - np.max(in_borders[1]) - 1
    y_min_border = np.min(in_borders[0])
    y_max_border = in_fits.shape[0] - np.max(in_borders[0]) - 1

    # Decide if it's worth removing the padding; if not, just return False output
    if ((x_min_border+x_max_border)<(0.1*in_fits.shape[1])) and ((y_min_border+y_max_border)<(0.1*in_fits.shape[0])):
        return band_dict

    # Slice smaller map out of original map
    out_fits = in_fits.copy()
    if y_min_border>0:
        out_fits = out_fits[y_min_border:,:]
    if x_min_border>0:
        out_fits = out_fits[:,x_min_border:]
    if y_max_border>0:
        out_fits = out_fits[:-y_max_border,:]
    if x_max_border>0:
        out_fits = out_fits[:,:-x_max_border]

    # Update header WCS to reflect changes
    out_wcs = astropy.wcs.WCS(naxis=2)
    out_wcs.wcs.crpix = [in_wcs.wcs.crpix[0]-x_min_border, in_wcs.wcs.crpix[1]-y_min_border]
    out_wcs.wcs.cdelt = in_wcs.wcs.cdelt
    out_wcs.wcs.crval = in_wcs.wcs.crval
    out_wcs.wcs.ctype = in_wcs.wcs.ctype
    out_header = out_wcs.to_header()

    # Save cutout to file
    out_fitspath = os.path.join( kwargs_dict['temp_dir_path'], 'Cutouts', source_dict['name'], source_dict['name']+'_'+band_dict['band_name']+'.fits' )
    out_cutout_hdu = astropy.io.fits.PrimaryHDU(data=out_fits, header=out_header)
    out_cutout_hdulist = astropy.io.fits.HDUList([out_cutout_hdu])
    out_cutout_hdulist.writeto(out_fitspath, clobber=True)

    # Repeat process for error map, if necessary
    if band_dict['use_error_map']==True:

        # Load in error map and extract WCS
        in_fits_error = astropy.io.fits.getdata(in_fitspath_error)

        # Slice smaller map out of original map
        out_fits_error = in_fits_error.copy()
        out_fits_error = out_fits_error[y_min_border:,:]
        out_fits_error = out_fits_error[:,x_min_border:]
        out_fits_error = out_fits_error[:-y_max_border,:]
        out_fits_error = out_fits_error[:,:-x_max_border]

        # Save cutout to file
        out_fitspath_error = os.path.join( kwargs_dict['temp_dir_path'], 'Cutouts', source_dict['name'], source_dict['name']+'_'+band_dict['band_name']+'_Error.fits' )
        out_cutout_hdu_error = astropy.io.fits.PrimaryHDU(data=out_fits_error, header=out_header)
        out_cutout_hdulist_error = astropy.io.fits.HDUList([out_cutout_hdu_error])
        out_cutout_hdulist_error.writeto(out_fitspath_error, clobber=True)

    # Return the directory of the newly-created cutout
    out_fitsdir = os.path.split(out_fitspath)[0]
    band_dict['band_dir'] = out_fitsdir
    return band_dict





# Define function that writes final aperture for given soruce to aperture output file
def RecordAperture(aperture_combined, source_dict, kwargs_dict):
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Writing apertures to file.'



    # Consturct string of line to be written
    aperture_string = str([ source_dict['name'], aperture_combined[0], aperture_combined[1], aperture_combined[2] ])#'name','semimaj_arcsec,axial_ratio,pos_angle\n'
    aperture_string = aperture_string.replace('[','').replace(']','').replace(' ','').replace('\'','')+'\n'

    # Write line to file
    aperture_table_file = open( kwargs_dict['aperture_table_path'], 'a')
    aperture_table_file.write(aperture_string)
    aperture_table_file.close()





# Define function that writes final aperture for given soruce to aperture output file
def RecordPhotom(photom_list, source_dict, bands_dict, kwargs_dict):
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Writing photometry to file.'



    # Find any bands not listed in the photom results, and add them as null measurements
    for band in bands_dict.keys():
        band_photom = False
        for photom_entry in photom_list:
            if band == photom_entry['band_name']:
                band_photom = True
                break
        if band_photom==True:
            continue
        elif band_photom==False:
            photom_null = {'band_name':band,
                            'ap_sum':np.NaN,
                            'ap_error':np.NaN}
            photom_list.append(photom_null)



    # Use band input table to establish order in which to put bands in results file, handling case of only one band
    photom_string = source_dict['name']
    bands_table = np.genfromtxt(kwargs_dict['bands_table_path'], delimiter=',', names=True, dtype=None)
    bands_list = bands_table['band_name']

    # Consturct string of line to be written
    if bands_list.shape==():
        bands_list = [bands_list.tolist()]
    for band_name in bands_list:
        for photom_entry in photom_list:
            if photom_entry['band_name']==band_name:
                photom_string += ','+str(photom_entry['ap_sum'])+','+str(photom_entry['ap_error'])
    photom_string += '\n'

    # Write line to file
    photom_table_file = open( kwargs_dict['photom_table_path'], 'a')
    photom_table_file.write(photom_string)
    photom_table_file.close()





# Define function that checks whether a decent amount of RAM is free before allowing things to progress
def MemCheck(pod, thresh_fraction=0.75, thresh_factor=20.0, swap_thresh_fraction=0.5, return_status=False):



    # Check whether psutil is available (as it is only for unix systems)
    try:
        import psutil
    except:
        if pod['verbose']:print '['+pod['id']+'] Library psutil not available (is this a Windows system?); unable to check RAM, proceeding regardless.'

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
            time.sleep(10.0+(5.0*np.random.rand()))
        else:
            break

        # Return memory status
        if return_status:
            return mem_stats





# Define function that makes a grid of aperture thumbnails for a given source's output
def ApertureThumbGrid(source_dict, bands_dict, kwargs_dict, aperture_list, aperture_combined):



    # Check that thumbnails are requested; if so, set up warnings and proceeed
    if kwargs_dict['thumbnails']!=True:
        return
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
    edge_long = np.floor(phi*float(edge_short))
    while int(edge_short*edge_long)<thumb_files:
        if int(edge_short+1)<int(edge_long):
            edge_short += 1
        else:
            edge_long += 1
    if (int(edge_short*edge_long)-int(edge_long))>=thumb_files:
        edge_short -= 1
    if (int(edge_short*edge_long)-int(edge_short))>=thumb_files:
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



    # Find largest beam size and outer annulus size, and hence work out thumbnail size that will contain the largest beam-convolved aperture
    beam_arcsec_max = 0.0
    outer_annulus_max = 0.0
    pix_arcsec_max = 0.0
    for band_name in bands_dict:
        if not os.path.exists(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'.fits')):
            continue
        band_pix_matrix = astropy.wcs.WCS(astropy.io.fits.getheader(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'.fits'))).pixel_scale_matrix
        band_pix_arcsec = 3600.0 * np.sqrt( np.min(np.abs(band_pix_matrix))**2.0 + np.max(np.abs(band_pix_matrix))**2.0 )
        if band_pix_arcsec>pix_arcsec_max:
            pix_arcsec_max = band_pix_arcsec
        if bands_dict[band_name]['beam_arcsec']>beam_arcsec_max:
            beam_arcsec_max = bands_dict[band_name]['beam_arcsec']
        if bands_dict[band_name]['annulus_outer']>outer_annulus_max:
            outer_annulus_max = bands_dict[band_name]['annulus_outer']
    thumb_rad = np.ceil( 1.0 * pix_arcsec_max ) + np.ceil( 1.2 * 0.5 * np.sqrt( (2.0*outer_annulus_max*aperture_combined[0])**2.0 + (beam_arcsec_max)**2.0 ) )
    pdb.set_trace()


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

        # Produce cutouts, and end loop
        if kwargs_dict['parallel']==True:
            thumb_pool.apply_async( ThumbCutout, args=(source_dict, bands_dict[band_name], kwargs_dict, img_input, thumb_rad,) )
        elif kwargs_dict['parallel']==False:
            ThumbCutout(source_dict, bands_dict[band_name], kwargs_dict, img_input, thumb_rad)
    thumb_pool.close()
    thumb_pool.join()



    # Begin main thumbnail plotting loop
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

        # Create and format image
        vars()['subfig'+str(b)] = aplpy.FITSFigure(img_output, figure=fig, subplot=[x_min, y_min, dx, dy])
        vars()['subfig'+str(b)].show_colorscale(cmap=bands_dict[band_name]['colour_map'], stretch='arcsinh', pmin=7.5, pmax=99.5)
        vars()['subfig'+str(b)].set_nan_color('black')
        vars()['subfig'+str(b)].axis_labels.hide()
        vars()['subfig'+str(b)].tick_labels.hide()
        vars()['subfig'+str(b)].ticks.hide()
        line_width = 4.0

        # Extract band-specific aperture dimensions
        band_ap_angle = aperture_list[b]['opt_angle']
        band_ap_axial_ratio = aperture_list[b]['opt_axial_ratio']
        band_ap_semimaj = (aperture_list[b]['opt_semimaj_arcsec'])/3600.0
        band_ap_semimin = band_ap_semimaj / band_ap_axial_ratio
        band_beam_width = bands_dict[band_name]['beam_arcsec'] / 3600.0
        band_ap_semimaj = ( band_ap_semimaj**2.0 + (0.5*band_beam_width)**2.0 )**0.5
        band_ap_semimin = ( band_ap_semimin**2.0 + (0.5*band_beam_width)**2.0 )**0.5
        band_ap_axial_ratio = band_ap_semimaj / band_ap_semimaj

        # Plot band-specific aperture (if one was provided)
        if isinstance(source_dict['aperture_bands_exclude'], str):
            aperture_bands_exclude = source_dict['aperture_bands_exclude'].split(';')
        else:
            aperture_bands_exclude = []
        if bands_dict[band_name]['consider_aperture']==True and band_name not in aperture_bands_exclude:
            vars()['subfig'+str(b)].show_ellipses(source_dict['ra'], source_dict['dec'], 2.0*band_ap_semimaj, 2.0*band_ap_semimin, angle=band_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width/2.0, linestyle='dotted')

        # Extract combined aperture dimensions
        comb_ap_angle = aperture_combined[2]
        comb_ap_axial_ratio = aperture_combined[1]
        comb_ap_semimaj = aperture_combined[0]/3600.0
        comb_ap_semimin = comb_ap_semimaj / comb_ap_axial_ratio
        comb_ap_semimaj = 0.5 * ( (2.0*comb_ap_semimaj)**2.0 + band_beam_width**2.0 )**0.5
        comb_ap_semimin = 0.5 * ( (2.0*comb_ap_semimin)**2.0 + band_beam_width**2.0 )**0.5

        # Plot combined aperture
        vars()['subfig'+str(b)].show_ellipses(source_dict['ra'], source_dict['dec'], 2.0*comb_ap_semimaj, 2.0*comb_ap_semimin, angle=comb_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width)

        # Plot combined background annulus
        band_ann_inner_semimaj = comb_ap_semimaj * 2.0 * bands_dict[band_name]['annulus_inner']
        band_ann_inner_semimin = comb_ap_semimin * 2.0 * bands_dict[band_name]['annulus_inner']
        band_ann_outer_semimaj = comb_ap_semimaj * 2.0 * bands_dict[band_name]['annulus_outer']
        band_ann_outer_semimin = comb_ap_semimin * 2.0 * bands_dict[band_name]['annulus_outer']
        vars()['subfig'+str(b)].show_ellipses(source_dict['ra'], source_dict['dec'], band_ann_inner_semimaj, band_ann_inner_semimin, angle=comb_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width/3.0, linestyle='dashed')
        vars()['subfig'+str(b)].show_ellipses(source_dict['ra'], source_dict['dec'], band_ann_outer_semimaj,band_ann_outer_semimin, angle=comb_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width/3.0, linestyle='dashed')

        # Plot label
        vars()['subfig'+str(b)].add_label(0.035, 0.92, bands_dict[band_name]['band_name'], relative=True, size=20, color='white', horizontalalignment='left')

        # Progress counters
        counter += 1
        column_counter += 1
        if np.mod(float(counter)+1, x_grid)==0:
            row_counter += 1
            column_counter = 0

    # Save figure, and remove temporary files
    fig.savefig( os.path.join(kwargs_dict['output_dir_path'],'Aperture_Fitting_Thumbnails',source_dict['name']+'_Thumbnail_Grid.png'), facecolor='white', dpi=100.0)
    [os.remove(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',processed_map)) for processed_map in os.listdir(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps')) if '.fits' in processed_map]
    fig.clear()
    plt.close('all')
    gc.collect()





# Define function that makes a grid of photometry thumbnails for a given source's output
def PhotomThumbGrid(source_dict, bands_dict, kwargs_dict):



    # Check that thumbnails are requested; if so, set up warnings and proceeed
    if kwargs_dict['thumbnails']!=True:
        return
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
        raise Exception('Aperture value caontains more than one entry for current galaxy')
    else:
        aperture_index = aperture_index[0][0]

    # Extract aperture corresponding to current source, dealing with special case where aperture file contains only one source
    if aperture_table['semimaj_arcsec'].shape==() and np.where( aperture_table['name']==source_dict['name'] )[0][0]==0:
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
        if int(edge_short+1)<int(edge_long):
            edge_short += 1
        else:
            edge_long += 1
    if (int(edge_short*edge_long)-int(edge_long))>=thumb_files:
        edge_short -= 1
    if (int(edge_short*edge_long)-int(edge_short))>=thumb_files:
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



    # Find largest beam size, outer annulus size, and pixel size - and hence work out thumbnail size that will contain the largest beam-convolved aperture
    beam_arcsec_max = 0.0
    outer_annulus_max = 0.0
    pix_arcsec_max = 0.0
    for band_name in bands_dict:
        if not os.path.exists(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'.fits')):
            continue
        band_pix_matrix = astropy.wcs.WCS(astropy.io.fits.getheader(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_name+'.fits'))).pixel_scale_matrix
        band_pix_arcsec = 3600.0 * np.sqrt( np.min(np.abs(band_pix_matrix))**2.0 + np.max(np.abs(band_pix_matrix))**2.0 )
        if band_pix_arcsec>pix_arcsec_max:
            pix_arcsec_max = band_pix_arcsec
        if bands_dict[band_name]['beam_arcsec']>beam_arcsec_max:
            beam_arcsec_max = bands_dict[band_name]['beam_arcsec']
        if bands_dict[band_name]['annulus_outer']>outer_annulus_max:
            outer_annulus_max = bands_dict[band_name]['annulus_outer']
    thumb_rad = np.ceil( 1.0 * pix_arcsec_max ) + np.ceil( 1.2 * 0.5 * np.sqrt( (2.0*outer_annulus_max*opt_semimaj_arcsec)**2.0 + (beam_arcsec_max)**2.0 ) )



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

        # Produce cutouts, and end loop
        if kwargs_dict['parallel']==True:
            thumb_pool.apply_async( ThumbCutout, args=(source_dict, bands_dict[band_name], kwargs_dict, img_input, thumb_rad,) )
        elif kwargs_dict['parallel']==False:
            ThumbCutout(source_dict, bands_dict[band_name], kwargs_dict, img_input, thumb_rad)
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
        comb_ap_semimaj = opt_semimaj_arcsec / 3600.0
        comb_ap_semimin = comb_ap_semimaj / comb_ap_axial_ratio
        comb_ap_semimaj = 0.5 * ( (2.0*comb_ap_semimaj)**2.0 + band_beam_width**2.0 )**0.5
        comb_ap_semimin = 0.5 * ( (2.0*comb_ap_semimin)**2.0 + band_beam_width**2.0 )**0.5

        # Plot combined aperture
        line_width = 4.0
        vars()['subfig'+str(w)].show_ellipses(source_dict['ra'], source_dict['dec'], 2.0*comb_ap_semimaj, 2.0*comb_ap_semimin, angle=comb_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width)

        # Plot combined background annulus
        band_ann_inner_semimaj = comb_ap_semimaj * 2.0 * bands_dict[band_name]['annulus_inner']
        band_ann_inner_semimin = comb_ap_semimin * 2.0 * bands_dict[band_name]['annulus_inner']
        band_ann_outer_semimaj = comb_ap_semimaj * 2.0 * bands_dict[band_name]['annulus_outer']
        band_ann_outer_semimin = comb_ap_semimin * 2.0 * bands_dict[band_name]['annulus_outer']
        vars()['subfig'+str(w)].show_ellipses(source_dict['ra'], source_dict['dec'], band_ann_inner_semimaj, band_ann_inner_semimin, angle=comb_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width/3.0, linestyle='dashed')
        vars()['subfig'+str(w)].show_ellipses(source_dict['ra'], source_dict['dec'], band_ann_outer_semimaj,band_ann_outer_semimin, angle=comb_ap_angle, edgecolor='#00FF40', facecolor='none', linewidth=line_width/3.0, linestyle='dashed')

        # Plot label
        vars()['subfig'+str(w)].add_label(0.035, 0.92, bands_dict[band_name]['band_name'], relative=True, size=20, color='white', horizontalalignment='left')

        # Progress counters
        counter += 1
        column_counter += 1
        if np.mod(float(counter)+1, x_grid)==0:
            row_counter += 1
            column_counter = 0

    # Save figure, and remove temporary files
    fig.savefig( os.path.join(kwargs_dict['output_dir_path'],'Photometry_Thumbnails',source_dict['name']+'_Thumbnail_Grid.png'), facecolor='white', dpi=100.0)
    [ os.remove(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',processed_map)) for processed_map in os.listdir(os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps')) if '.fits' in processed_map ]
    fig.clear()
    plt.close('all')
    gc.collect()





# Define function for producing an appropriately-sized cutout in current band
def ThumbCutout(source_dict, band_dict, kwargs_dict, img_input, thumb_rad):



    # Construct output path
    img_output = os.path.join(kwargs_dict['temp_dir_path'],'Processed_Maps',source_dict['name']+'_'+band_dict['band_name']+'_Thumbnail.fits')

    # Produce cutout, and neaten header as necessary (to quell noisy APLpy verbosity later on)
    cutout_parallel = ( not kwargs_dict['parallel'] )
    cutout_hdulist = ChrisFuncs.FitsCutout( img_input, source_dict['ra'], source_dict['dec'], thumb_rad, exten=0, variable=True, reproj=True, parallel=cutout_parallel, fast=True )
    cutout_data = cutout_hdulist[0].data
    cutout_header = cutout_hdulist[0].header
    cutout_header['EQUINOX'] = 2000.0
    cutout_header['EPOCH'] = 2000.0

    # Write thumbnail cutout to file, and end output supression
    astropy.io.fits.writeto(img_output, cutout_data, header=cutout_header, clobber=True)
