# Import smorgasbord
import os
import sys
sys.path.insert(0, '../')
import gc
import pdb
import time
import csv
import psutil
import numpy as np
import astropy.io.fits
import astropy.wcs
import ChrisFuncs



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

    # make sure appropriate cutout sub-directories exist in temp directory
    if not os.path.exists( os.path.join( temp_dir_path, 'Cutouts' ) ):
        os.mkdir( os.path.join( temp_dir_path, 'Cutouts' ) )
    if not os.path.exists( os.path.join( temp_dir_path, 'Cutouts', source_dict['name'] ) ):
        os.mkdir( os.path.join( temp_dir_path, 'Cutouts', source_dict['name'] ) )

    # Using standard filename format, construct full file path, and work out whether the file extension is .fits or .fits.gz
    in_fitspath = os.path.join( band_dict['band_dir'], source_dict['name']+'_'+band_dict['band_name'] )
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
def MemCheck(pod, thresh_fraction=0.8, thresh_factor=10.0, return_status=False):

    # Start infinite loop
    wait_initial = True
    wait_count = 0
    while True:
        mem_wait = False

        # Assess how much RAM is free
        mem_stats = psutil.virtual_memory()
        mem_usage = 1.0-(float(mem_stats[1])/float(mem_stats[0]))

        # Require wait if less than 20% of RAM is free
        if thresh_fraction!=None:
            if mem_usage>=float(thresh_fraction):
                mem_wait = True

        # Also, require wait if the amount of RAM free is more than some multiple the size of the current file
        mem_free = float(psutil.virtual_memory()[4])
        if thresh_factor!=None:
            if ( float(thresh_factor) * float(pod['in_fitspath_size']) )>mem_free:
                if pod['in_fitspath_size']!=None:
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
            time.sleep(5.0+(10.0*np.random.rand()))
        else:
            break

        # Return memory status
        if return_status:
            return mem_stats

