# Import the relevant PTS classes and modules
import sys
import os
sys.path.insert(0, '../')
sys.path.append( str( os.path.join( os.path.split( os.path.dirname(os.path.abspath(__file__)) )[0], 'CAAPR_AstroMagic', 'PTS') ) )
import numpy as np
import warnings
import signal
import shutil
import pts
import pdb
import pickle
import time as pytime
import astropy.io.fits
import astropy.wcs
from astropy.units import Unit
import pyregion
import ChrisFuncs
from pts.magic.core.image import Image
from pts.magic.misc.imageimporter import ImageImporter
from pts.magic.sources.finder import SourceFinder
from pts.magic.catalog.importer import CatalogImporter
from pts.magic.sources.extractor import SourceExtractor
from pts.core.tools import configuration
from pts.core.tools import logging, time, filesystem
from pts.magic.basics.region import Region





# Define function that wraps the AstroMagic calling procedure
def Magic(pod, source_dict, band_dict, kwargs_dict):
    in_fitspath = pod['in_fitspath']
    temp_dir_path = pod['temp_dir_path']

    # Ensure star-subtraction is actually required
    if kwargs_dict['starsub']!=True:
        return pod
    if band_dict['remove_stars']!=True:
        return pod
    if band_dict['band_name'] in str(source_dict['starsub_bands_exclude']).split(';'):
        print '['+pod['id']+'] User explicitly excluded current band from star subtraction for this source.'
        return pod
    if str(source_dict['starsub_bands_exclude'])=='True':
        print '['+pod['id']+'] User explicitly requested no star subtraction for this source.'
        return pod
    if pod['verbose']: print '['+pod['id']+'] Removing foreground stars and background galaxies with PTS AstroMagic.'



    # The paths to the image and output (absolute or relative to the current working directory)
    output_path = os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name'])
    image_path = in_fitspath



    # Setting the log level to "ERROR" disables all output from PTS except error messages (probably want to see those); full options are "DEBUG", "INFO", "WARNING", "ERROR", "SUCCESS"
    logging.setup_log(level="ERROR")



    # The path to the bad region (as explained in the previous script that I sent you)
    bad_region_path = None

    # The FWHM of the image (if known)
    fwhm = band_dict['beam_arcsec'] * Unit("arcsec")

    # Import the image
    importer = ImageImporter()
    importer.run(image_path, bad_region_path=bad_region_path, fwhm=fwhm, find_error_frame=False)

    # Get the imported image
    image = importer.image

    # Get the mask of bad pixels
    bad_mask = image.masks.bad if "bad" in image.masks else None



    # If version of cutout alrady processed by AstroMagic is present, use it; else, commence regular processing
    if os.path.exists( os.path.join( temp_dir_path, 'AstroMagic', band_dict['band_name'], source_dict['name']+'_'+band_dict['band_name']+'_StarSub.fits') ):
        if pod['verbose']: print '['+pod['id']+'] AstroMagic accessing pre-processed data for this map.'
        am_output = astropy.io.fits.getdata( os.path.join( temp_dir_path, 'AstroMagic', band_dict['band_name'], source_dict['name']+'_'+band_dict['band_name']+'_StarSub.fits') )
        pod['cutout'] = am_output
        pod['pre_reg'] = True
        return pod
    else:

        # The path to the directory where all the output will be placed
        if os.path.exists(os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name'])):
            shutil.rmtree(os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name']))
            os.mkdir(os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name']))
        else:
            os.mkdir(os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name']))

        # Make copies of pre-fetched catalogues, to prevent simultaneous access conflicts
        shutil.copy(os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic', 'Stars.cat'), os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name'], 'Stars.cat'))
        shutil.copy(os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic', 'Galaxies.cat'), os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name'], 'Galaxies.cat'))



        # Create a CatalogImporter instance, and run it to import catalogues
        catalog_importer = CatalogImporter()
        catalog_importer.config.stars.use_catalog_file = True
        catalog_importer.config.galaxies.use_catalog_file = True
        catalog_importer.config.stars.catalog_path = os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name'], 'Stars.cat')
        catalog_importer.config.galaxies.catalog_path = os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name'], 'Galaxies.cat')
        catalog_importer.run(image.frames.primary)

        # If currently handling cut-down thumbnail, only use relevant portion of catalogues
        if 'starsub_thumbnail' in pod.keys():
            if pod['starsub_thumbnail']==True:
                galaxy_catalog_thumb = ThumbCatalogue(pod, source_dict, band_dict, kwargs_dict, catalog_importer)
                catalog_importer.galactic_catalog = galaxy_catalog_thumb



        # Create a SourceFinder instance
        finder = SourceFinder()

        # If you don't want to do the 'find_other_sources' step, comment-out the line below
        finder.config.find_other_sources = False # default is True

        # Downsample map for faster run time
        if int(band_dict['downsample_factor'])>1:
            if pod['verbose']: print '['+pod['id']+'] AstroMagic will run on copy of map downsampled by factor of '+str(int(band_dict['downsample_factor']))+', to improve speed.'
            finder.config.downsample_factor = int( band_dict['downsample_factor'] )

        # Run the source finder
        if pod['verbose']: print '['+pod['id']+'] AstroMagic locating online catalogue sources in map.'
        special_region = None # Not important except for debugging
        ignore_region = None # Not important except when certain areas need to be completely ignored from the extraction procedure
        finder.config.build_catalogs = False # For using pre-fetched catalogue files
        try:
            finder.run(image.frames.primary, catalog_importer.galactic_catalog, catalog_importer.stellar_catalog, special_region, ignore_region, bad_mask)
        except RuntimeError as e:
            if 'found' in e.message:
                if pod['verbose']: print '['+pod['id']+'] AstroMagic found no sources to remove.'
                pod['pre_reg'] = True
                return pod
            else:
                pdb.set_trace()



        # Save the galaxy region
        galaxy_region = finder.galaxy_region
        galaxy_region_path = filesystem.join(output_path, "galaxies.reg")
        galaxy_region.save(galaxy_region_path)

        # Save the star region
        star_region = finder.star_region
        star_region_path = filesystem.join(output_path, "stars.reg")
        if star_region is not None: star_region.save(star_region_path)

        # Save the saturation region
        saturation_region = finder.saturation_region
        saturation_region_path = filesystem.join(output_path, "saturation.reg")
        if saturation_region is not None: saturation_region.save(saturation_region_path)

        # Save the region of other sources
        other_region = finder.other_region
        other_region_path = filesystem.join(output_path, "other_sources.reg")
        if other_region is not None: other_region.save(other_region_path)



        # Get the segmentation maps (galaxies, stars and other sources) from the SourceFinder
        galaxy_segments = finder.galaxy_segments
        star_segments = finder.star_segments
        other_segments = finder.other_segments

        # Make sure target galaxy isn't identified as star segment
        target_segment = star_segments[ pod['centre_i'], pod['centre_j'] ]
        star_segments[ np.where( star_segments==target_segment ) ] = 0.0




        # Handle stars that have been conflated with the target galaxy
        star_segments = OverlargeStars(pod, star_segments, saturation_region_path, star_region_path, galaxy_region_path, image, source_dict, band_dict, temp_dir_path)

        # Region files can be adjusted by the user; if this is done, they have to be reloaded
        star_region = Region.from_file(star_region_path.replace('.reg','_revised.reg'))
        saturation_region = Region.from_file(saturation_region_path.replace('.reg','_revised.reg'))
        """
        # Remove all but target galaxy from galaxy region file
        ExcessGalaxies(galaxy_region_path, galaxy_principal)
        galaxy_region = Region.from_file(galaxy_region_path.replace('.reg','_revised.reg'))
        """
        # Remopve all galaxies from galaxy region file
        shutil.copy2(galaxy_region_path, galaxy_region_path.replace('.reg','_revised.reg'))
        gal_header = '# Region file format: DS9 version 4.1\n'
        gal_file_new = open( galaxy_region_path.replace('.reg','_revised.reg'), 'w')
        gal_file_new.write(gal_header)
        gal_file_new.close()
        galaxy_region = Region.from_file(galaxy_region_path.replace('.reg','_revised.reg'))



        # Create a map of the the segmentation maps
        segments = Image("segments")

        # Add the segmentation map of the galaxies
        segments.add_frame(galaxy_segments, "galaxies")

        # Add the segmentation map of the saturated stars
        if star_segments is not None: segments.add_frame(star_segments, "stars")

        # Add the segmentation map of the other sources
        if other_segments is not None: segments.add_frame(other_segments, "other_sources")

        # Save the FITS file with the segmentation maps
        path = filesystem.join(output_path, "segments.fits")
        segments.save(path)



        # Create an SourceExtractor instance
        if pod['verbose']: print '['+pod['id']+'] AstroMagic extracting background sources.'
        extractor = SourceExtractor()

        # Run the source extractor
        extractor.run(image.frames.primary, galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments)



        # Determine the path to the result
        result_path = filesystem.join( temp_dir_path, 'AstroMagic', band_dict['band_name'], source_dict['name']+'_'+band_dict['band_name']+'_StarSub.fits' )

        # Save the resulting image as a FITS file
        image.frames.primary.save(result_path, header=image.original_header)



        # Grab AstroMagic output and return in pod
        am_output = astropy.io.fits.getdata( result_path )
        pod['cutout'] = am_output
        pod['pre_reg'] = True
        return pod





# Define function to saturation regions that wholly intersect target galaxy
def OverlargeStars(pod, star_segments, sat_path, star_path, gal_path, image, source_dict, band_dict, temp_dir_path):



    # Check if there are any saturation regions in saturation file; if not, return un-altered saturation files
    sat_line_count = 0
    with open(sat_path) as sat_file:
        sat_line_count = sum(1 for _ in sat_file)
    if sat_line_count<2:
        shutil.copy2(sat_path, sat_path.replace('.reg','_revised.reg'))
        shutil.copy2(star_path, star_path.replace('.reg','_revised.reg'))
        return



    # Read in saturation region file, and loop over each entry, recording area and index
    try:
        sat_regions = pyregion.open(sat_path)
        sat_indices, sat_areas = [], []
        for sat_region in sat_regions:
            sat_indices.append( float(sat_region.attr[1]['text']) )
            sat_areas.append( np.pi * float(sat_region.coord_list[2]) * float(sat_region.coord_list[3]) )
        #sat_array = np.array([sat_areas, sat_indices]).transpose()
        #sat_indices_out = sat_array[:,1].astype(int).tolist()
    except ValueError as error_message:
        if error_message.message=='need more than 0 values to unpack':
            sat_regions = []
    except:
        pdb.set_trace()

    # Open galaxy catalogue file, and determine the "primary" name of the one that has been deemed principal
    gal_cat = astropy.table.Table.read(os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name'], 'Galaxies.cat'), format='ascii')
    where_principal = np.where( np.array([ gal_cat['Name'][i].replace(' ','')==source_dict['name'] for i in range(0,len(gal_cat)) ])==True )
    if where_principal[0].shape[0]==0:
        gal_principal = 'NULL'
    else:
        gal_principal = gal_cat['Name'][where_principal][0].replace(' ','')

    # Loop over galaxy region file, identify region corresponding to target galaxy, and create mask array of it
    gal_regions = pyregion.open(gal_path)
    gal_found = False
    for gal_region in gal_regions:
        if 'text' in gal_region.attr[1].keys():
            gal_name = gal_region.attr[1]['text'].replace(' ','').replace(' (principal)','')
            if gal_name==gal_principal:
                gal_found = True
                break
    if gal_found==False:
        shutil.copy2(sat_path, sat_path.replace('.reg','_revised.reg'))
        shutil.copy2(star_path, star_path.replace('.reg','_revised.reg'))
        return star_segments



    # Loop back over the saturation regions, not keeping those that aren't being retained, or which wholly encompass target galaxy, or are too lose to galaxy centre
    sat_regions_out = pyregion.ShapeList([])
    principal_dist_beam_thresh = 2.0
    star_regions = pyregion.open(star_path)
    for sat_region in sat_regions:

        # Do initial check if saturation region encompasses centre of principal galaxy
        i_dist = gal_region.coord_list[0] - sat_region.coord_list[0]
        j_dist = gal_region.coord_list[1] - sat_region.coord_list[1]
        if (i_dist**2.0+j_dist**2.0)<=(np.min(sat_region.coord_list[2:4]))**2.0:

            # If initial check indicates risk, do full mask check
            sat_mask = ChrisFuncs.Photom.EllipseMask( np.zeros(image.shape), np.max(sat_region.coord_list[2:4]), np.max(sat_region.coord_list[2:4])/np.min(sat_region.coord_list[2:4]), sat_region.coord_list[4], sat_region.coord_list[1], sat_region.coord_list[0])
            if sat_mask[ pod['centre_i'], pod['centre_j'] ]==1.0:

                # If saturation region envelops galaxy core, find corresponding star region
                sat_region_index = sat_region.attr[1]['text']
                for star_region in star_regions:
                    if 'text' in star_region.attr[1].keys():
                        if star_region.attr[1]['text']==sat_region_index:

                            # Use star region info to create more sensible saturationr egion
                            star_i_dist = gal_region.coord_list[0] - star_region.coord_list[0]
                            star_j_dist = gal_region.coord_list[1] - star_region.coord_list[1]
                            sat_dist_new = np.max([ 0, np.sqrt(star_i_dist**2.0+star_j_dist**2.0)-np.max(gal_region.coord_list[2:4]) ])
                            sat_region.coord_list[0] = star_region.coord_list[0]
                            sat_region.coord_list[1] = star_region.coord_list[1]
                            sat_region.coord_list[2] = sat_dist_new
                            sat_region.coord_list[3] = sat_dist_new

                            # Update saturated star segments map
                            star_segments[ np.where(star_segments==float(sat_region_index)) ] = 0.0
                            sat_mask_new = ChrisFuncs.Photom.EllipseMask( np.zeros(image.shape), np.max(sat_region.coord_list[2:4]), np.max(sat_region.coord_list[2:4])/np.min(sat_region.coord_list[2:4]), sat_region.coord_list[4], sat_region.coord_list[1], sat_region.coord_list[0])
                            star_segments[ np.where(sat_mask_new==1.0) ] = float(sat_region_index)
                            sat_regions_out.append(sat_region)
                            #sat_indices_bad.append( float(sat_region.attr[1]['text']) )
                            continue
        else:
            sat_regions_out.append(sat_region)



    # Now read in and loop over the star region file
    star_regions = pyregion.open(star_path)
    star_regions_out = pyregion.ShapeList([])
    for star_region in star_regions:

        # Remove stars that are located too close to centre of principal galaxy
        principal_dist = np.sqrt( (gal_region.coord_list[0]-star_region.coord_list[0])**2.0 + (gal_region.coord_list[1]-star_region.coord_list[1])**2.0 )
        if principal_dist<=((band_dict['beam_arcsec']/pod['pix_arcsec'])*principal_dist_beam_thresh):
            continue

        # If star has passed all criteria, record it to output
        star_regions_out.append(star_region)



    # Save updated regions to file
    if len(sat_regions_out)>0:
        sat_regions_out.write(sat_path.replace('.reg','_revised.reg'))
    else:
        sat_regions_out = open(sat_path.replace('.reg','_revised.reg'), 'w')
        sat_regions_out.close()
    if len(star_regions_out)>0:
        star_regions_out.write(star_path.replace('.reg','_revised.reg'))
    else:
        star_regions_out = open(star_path.replace('.reg','_revised.reg'), 'w')
        star_regions_out.close()

    # Return updated segments map
    if gal_principal=='NULL':
        gal_principal = None
    return star_segments





# Define function that removed all except target galaxy from galaxy region file
def ExcessGalaxies(gal_region_path, gal_principal):



    # Read in galaxy region file, and loop over each entry
    gal_principal_found = False
    if gal_principal!=None:
        with open(gal_region_path) as gal_file:
            for gal_line in gal_file:

                # Skip over point markers; only interested in galaxy ellipses
                if 'image;ellipse' not in gal_line:
                    continue

                # Identify if current line refers to target galaxy
                if '{'+gal_principal+'}' in gal_line.replace(' ',''):
                    gal_principal_found = True

                    # If target galaxy region file, regord line
                    gal_region = gal_line

                    # Also create point marker for galaxy
                    gal_point = gal_line
                    gal_point = gal_point.split(',')[0]+','+gal_point.split(',')[1]+') # point = x\n'
                    gal_point = gal_point.replace('image;ellipse','image;point')

                    # Create new copy of galaxy region file, containing only target galaxy
                    shutil.copy2(gal_region_path, gal_region_path.replace('.reg','_revised.reg'))
                    gal_header = '# Region file format: DS9 version 4.1\n'
                    gal_file_new = open( gal_region_path.replace('.reg','_revised.reg'), 'w')
                    gal_file_new.write(gal_header)
                    gal_file_new.write(gal_region)
                    gal_file_new.write(gal_point)
                    gal_file_new.close()

    # Recrd if principal galaxy not found
    if gal_principal_found==False:
        gal_principal = None

    # If principal not found, of no principal galaxy name provided, just produce empty galaxy region file
    if gal_principal==None:
        shutil.copy2(gal_region_path, gal_region_path.replace('.reg','_revised.reg'))
        gal_header = '# Region file format: DS9 version 4.1\n'
        gal_file_new = open( gal_region_path.replace('.reg','_revised.reg'), 'w')
        gal_file_new.write(gal_header)
        gal_file_new.close()





# Define function to only deal with relevant portion of imported catalogues, when only working with thumbnail cutouts
def ThumbCatalogue(pod, source_dict, band_dict, kwargs_dict, catalog_importer):



    # Make copy of galaxy catalogue, and make arrays to record thumb coverage
    gal_cat = catalog_importer.galactic_catalog.copy()
    gal_thumb = np.array( [False]*gal_cat.as_array().shape[0] )
    gal_centre = np.array( pod['cutout'].shape )/2.0
    gal_centre_dist_min = 1E50

    # Loop over each entry in galaxy catalogue
    for i in range(0, gal_cat.as_array().shape[0]):
        gal_xy = pod['in_wcs'].wcs_world2pix( np.array([[ gal_cat['Right ascension'][i], gal_cat['Declination'][i] ]]), 0 )

        # Check if coord lies outside thumbnail
        if gal_xy[0][1]<=0:
            continue
        if gal_xy[0][1]>=pod['cutout'].shape[0]:
            continue
        if gal_xy[0][0]<=0:
            continue
        if gal_xy[0][0]>=pod['cutout'].shape[1]:
            continue

        # Take note if coord lies within thumbnail
        gal_thumb[i] = True

        # Record if current galaxy is closest yet to centre
        gal_centre_dist = np.sqrt( (gal_xy[0][1]-float(pod['cutout'].shape[0]))**2.0 + (gal_xy[0][0]-float(pod['cutout'].shape[1]))**2.0 )
        if gal_centre_dist<gal_centre_dist_min:
            gal_centre_dist_min = gal_centre_dist
            gal_centre = np.array( [False]*gal_cat.as_array().shape[0] )
            gal_centre[i] = True

    # Reduce catalogue to only those entries that lie within thumbnail cutout; if this removes all galaxies, make dummy entry and continue
    gal_cat = gal_cat[ np.where( gal_thumb==True ) ]
    if gal_cat.as_array().shape[0]==0:
        gal_cat.add_row()
        for col in gal_cat.colnames:
            gal_cat[col].mask[0] = True
        gal_cat[0]['Name'] = source_dict['name']
        gal_cat[0]['Right ascension'] = source_dict['ra']
        gal_cat[0]['Declination'] = source_dict['dec']
        gal_cat[0]['Principal'] = True
        return gal_cat

    # If no galaxies are now labelled as principal, set the most central galaxy to be
    if np.where( gal_cat['Principal']==True )[0].shape[0]==0:
        gal_centre = gal_centre[ np.where( gal_thumb==True ) ]
        gal_cat['Principal'] = gal_centre

    # Return  new catalogue
    return gal_cat





# Define function that acquires online catalogues for astromagic processing
def PreCatalogue(source_dict, bands_dict, kwargs_dict):



    # If absoltuely no star subtractionis required for this source, immediately return
    if kwargs_dict['starsub']==False:
        return
    if str(source_dict['starsub_bands_exclude'])=='True':
        return

    # If star subtraction is possibly required, check band-by-band
    if kwargs_dict['starsub']==True:
        star_sub_check = False
        for band in bands_dict.keys():
            if bands_dict[band]==None:
                continue
            if bands_dict[band]['remove_stars']==True:
                star_sub_check = True
        if star_sub_check==False:
            return

    # If all checks passed, and star subtraction is required, inform user and make sure that AstroMagic temp directory is clear
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] PTS AstroMagic retrieving list of foreground stars in map from online catalogues.'
    if os.path.exists(os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic')):
        shutil.rmtree(os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic'))
    os.mkdir(os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic'))



    # Loop over each band, determining map sizes
    diam_max = 0.0
    for band in bands_dict.keys():
        band_dict = bands_dict[band]

        # Check that this band requires star subtraction
        if band_dict['remove_stars']!=True:
            continue

        # Determine fits path
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
        if file_found==False:
            continue

        # Work out map size
        band_header = astropy.io.fits.getheader(in_fitspath)
        band_wcs = astropy.wcs.WCS(band_header)
        band_cdelt = band_wcs.wcs.cdelt.max()
        diam = np.max([ band_cdelt*float(band_header['NAXIS1']), band_cdelt*float(band_header['NAXIS2']) ])
        if diam>diam_max:
            file_max = in_fitspath
            diam_max = diam



    # Run catalogue-prefetching in a try statement, to catch a pernicious error
    try_counter = 0
    try_success = False
    while try_counter<10:
        if try_success==True:
            break
        try:

            # Get AstroMagic catalogue object reference fits
            logging.setup_log(level="ERROR")
            importer = ImageImporter()
            importer.run(file_max, find_error_frame=False)
            image = importer.image

            # Run catalogue importer on dummy fits, and save results
            catalog_importer = CatalogImporter()
            catalog_importer.config.writing.galactic_catalog_path = os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic', 'Galaxies.cat')
            catalog_importer.config.writing.stellar_catalog_path = os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic', 'Stars.cat')
            catalog_importer.run(image.frames.primary)
            catalog_importer.write_galactic_catalog()
            catalog_importer.write_stellar_catalog()

            # Record success and break
            try_success = True
            break

        # Handle error
        except ValueError as e:
            if 'rebin' in e.message:
                if kwargs_dict['verbose']: print '['+source_dict['name']+'] PTS AstroMagic encountered an error whilst pre-fetching stellar catalogues; re-attemping.'
                try_counter += 1
            else:
                pdb.set_trace()
    if try_counter>=10:
        pdb.set_trace()






"""
# Use IRSA header template service to create dummy header
try:
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] Accessing IRSA header template service to create generic template header.'
    dummy_cdelt_arcsec = (diam_max/dummy_pix)*3600.0
    url = 'http://irsa.ipac.caltech.edu/cgi-bin/HdrTemplate/nph-hdr?location='+str(source_dict['ra'])+'%2C+'+str(source_dict['dec'])+'&system=Equatorial&equinox=2000.&width='+str(diam_max)+'&height='+str(diam_max)+'&resolution='+str(dummy_cdelt_arcsec)+'&rotation=0.0'
    sys.stdout = open(os.devnull, "w")
    ChrisFuncs.wgetURL(url, os.path.join(kwargs_dict['temp_dir_path'],'Header_Template.txt'), clobber=True, auto_retry=False)
    sys.stdout = sys.__stdout__

# Produce dummy header
except:
dummy_cdelt = diam_max / dummy_pix
#dummy_cdelt_arcsec = (diam_max/100.0)*3600.0
dummy_header_in = open(os.path.join( os.path.split( os.path.dirname(os.path.abspath(__file__)) )[0], 'CAAPR_AstroMagic','Header_Template.txt'), 'r')
dummy_header_string = dummy_header_in.read()
dummy_header_in.close()
dummy_header_string = dummy_header_string.replace('RA_PLACEHOLDER', str(source_dict['ra']))
dummy_header_string = dummy_header_string.replace('DEC_PLACEHOLDER', str(source_dict['dec']))
dummy_header_string = dummy_header_string.replace('CDELT1_PLACEHOLDER', str(dummy_cdelt))
dummy_header_string = dummy_header_string.replace('CDELT2_PLACEHOLDER', str(dummy_cdelt))
dummy_header_out = open( os.path.join( kwargs_dict['temp_dir_path'],'Header_Template.txt' ), 'w')
dummy_header_out.write(dummy_header_string)
dummy_header_out.close()

# Create dummy map
dummy_header = astropy.io.fits.Header.fromfile( os.path.join(kwargs_dict['temp_dir_path'],'Header_Template.txt'), sep='\n', endcard=False, padding=False)
dummy_map = np.zeros([ int(dummy_header['NAXIS1']), int(dummy_header['NAXIS2']) ])
#dummy_cdelt_arcsec = (diam_max/float(dummy_header['NAXIS2']))*3600.0
dummy_file = os.path.join( kwargs_dict['temp_dir_path'],'FITS_Template.fits' )
astropy.io.fits.writeto( dummy_file, dummy_map, header=dummy_header)
"""

"""
# Save catalogues to file
catalog_importer.galactic_catalog.write(os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic', 'Galaxies.cat'), format='ascii.commented_header')
catalog_importer.stellar_catalog.write(os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic', 'Stars.cat'), format='ascii.commented_header')

# Run AstroMagic finder on dummy
finder = SourceFinder()
finder.config.writing.galactic_catalog_path = os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic', 'Galaxies.cat')
finder.config.writing.stellar_catalog_path = os.path.join(kwargs_dict['temp_dir_path'], 'AstroMagic', 'Stars.cat')
finder.config.find_other_sources = False
special_region = None
ignore_region = None
bad_mask = image.masks.bad if 'bad' in image.masks else None
finder.run(image.frames.primary, catalog_importer.galactic_catalog, catalog_importer.stellar_catalog, special_region, ignore_region, bad_mask)

# Export catalogue files
finder.write_galactic_catalog()
finder.write_stellar_catalog()

# Return catalogue importer object
return catalog_importer
"""

"""
# Define how stars not fit by PSF should be masked
finder.config.stars.source_psf_sigma_level = 3.0 # 4.0 is the default, change this to whatever value you want
finder.config.stars.fwhm.scale_factor = 2.0 # 1.0 is the default, change this to 2.0 for example
finder.config.stars.fwhm.measure = 'max' # Can change to median, mean, or max
finder.config.stars.saturation.dilate = True
finder.config.stars.saturation.iterations = 1
finder.config.stars.saturation.connectivity = 1
"""