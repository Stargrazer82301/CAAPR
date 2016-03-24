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
def Magic(pod, source_dict):
    pod_in = pod
    in_fitspath = pod['in_fitspath']
    temp_dir_path = pod['temp_dir_path']
    band_dict = pod['band_dict']
    source_dict = pod['source_dict']



    # Handle repeated attempts and timeouts
    def Handler(signum, frame):
        raise Exception("Timout!")
    signal.signal(signal.SIGALRM, Handler)
    try:
        panic
        signal.alarm(3600)



        # The path to the image (absolute or relative to the current working directory)
        image_path = in_fitspath

        # The path to the directory where all the output will be placed
        if os.path.exists(os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name'])):
            shutil.rmtree(os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name']))
        os.mkdir(os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name']))
        output_path = os.path.join(temp_dir_path, 'AstroMagic', band_dict['band_name'])



        # Setting the log level to "ERROR" disables all output from PTS except error messages (probably want to see those); full options are "DEBUG", "INFO", "WARNING", "ERROR", "SUCCESS"
        logging.setup_log(level="ERROR")



        # The path to the bad region (as explained in the previous script that I sent you)
        bad_region_path = None

        # The FWHM of the image (if known)
        fwhm = band_dict['beam_arcsec'] * Unit("arcsec")

        # Import the image
        importer = ImageImporter()
        importer.run(image_path, bad_region_path=bad_region_path, fwhm=fwhm)

        # Get the imported image
        image = importer.image

        # Get the mask of bad pixels
        bad_mask = image.masks.bad if "bad" in image.masks else None



    #        # Use pre-imported catalogue object
    #        catalog_importer = source_dict['pre_catalogue']

        # Create a CatalogImporter instance
        if pod['verbose']: print '['+pod['id']+'] AstroMagic retrieving list of sources in map from online catalogues.'
        catalog_importer = CatalogImporter()

        # Run the catalog importer
        catalog_importer.run(image.frames.primary)



        # Create a SourceFinder instance
        finder = SourceFinder()

        # If you don't want to do the "find_other_sources" step, comment-out the line below
        finder.config.find_other_sources = False # default is True



        # Run the source finder
        if pod['verbose']: print '['+pod['id']+'] AstroMagic locating online catalogue sources in map.'
        special_region = None # Not important except for debugging
        ignore_region = None # Not important except when certain areas need to be completely ignored from the extraction procedure
        finder.run(image.frames.primary, catalog_importer.galactic_catalog, catalog_importer.stellar_catalog, special_region, ignore_region, bad_mask)



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



        # Create an image with the segmentation maps
        if pod['verbose']: print '['+pod['id']+'] AstroMagic generating segmentation maps.'
        segments = Image("segments")

        # Get the segmentation maps (galaxies, stars and other sources) from the SourceFinder
        galaxy_segments = finder.galaxy_segments
        star_segments = finder.star_segments
        other_segments = finder.other_segments

        # Add the segmentation map of the galaxies
        segments.add_frame(galaxy_segments, "galaxies")

        # Add the segmentation map of the saturated stars
        if star_segments is not None: segments.add_frame(star_segments, "stars")

        # Add the segmentation map of the other sources
        if other_segments is not None: segments.add_frame(other_segments, "other_sources")

        # Save the FITS file with the segmentation maps
        path = filesystem.join(output_path, "segments.fits")
        segments.save(path)



        # Only process the most conspicuous foreground stars, to save time
        BrightestStars(saturation_region_path, star_region_path, galaxy_region_path, image, source_dict, percentile=80.0, maxtot=50)

        # Region files can be adjusted by the user; if this is done, they have to be reloaded
        star_region = Region.from_file(star_region_path.replace('.reg','_revised.reg'))
        saturation_region = Region.from_file(saturation_region_path.replace('.reg','_revised.reg'))
        if finder.config.find_other_sources==True: other_region = Region.from_file(other_region_path)



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
        am_output = astropy.io.fits.getdata( os.path.join( temp_dir_path, 'AstroMagic', band_dict['band_name'], source_dict['name']+'_'+band_dict['band_name']+'_StarSub.fits') )
        pod['cutout'] = am_output
        signal.alarm(0)
        return pod



    # Handle failure
    except:
        signal.alarm(0)
        print '['+pod['id']+'] AstroMagic processing failed; retaining standard image.'
        pod = pod_in
        return pod





# Define function to select only the most prominent foreground stars identified by AstroMagi
def BrightestStars(sat_path, star_path, gal_path, image, source_dict, percentile=80.0, maxtot=50):



    # Read in saturation region file, and loop over each entry, recording area and index
    sat_regions = pyregion.open(sat_path)
    sat_indices, sat_areas = [], []
    for sat_region in sat_regions:
        sat_indices.append( float(sat_region.attr[1]['text']) )
        sat_areas.append( np.pi * float(sat_region.coord_list[2]) * float(sat_region.coord_list[3]) )

    # Determine the area that corresponds to the requested percentile cut (or max total), and keep only saturation regions that satisfy those criteria
    sat_array = np.array([sat_areas, sat_indices]).transpose()
    sat_cutoff = np.percentile( sat_array[:,0], float(percentile) )
    if np.where(sat_areas>sat_cutoff)[0].shape[0]>maxtot:
        sat_cutoff = np.sort(sat_array)[::-1][int(maxtot)-1]
    sat_array = sat_array[ np.where( sat_array[:,0]>=sat_cutoff ) ]
    sat_indices_out = sat_array[:,1].astype(int).tolist()

    # Loop over galaxy region file, identify region corresponding to target galaxy, and create mask array of it
    gal_regions = pyregion.open(gal_path)
    for gal_region in gal_regions:
        if 'text' in gal_region.attr[1].keys():
            gal_name = gal_region.attr[1]['text'].replace(' ','')
            if gal_name==source_dict['name']:
                break
    gal_mask = ChrisFuncs.Photom.EllipseMask( np.zeros(image.shape), np.max(gal_region.coord_list[2:4]), np.max(gal_region.coord_list[2:4])/np.min(gal_region.coord_list[2:4]), gal_region.coord_list[4], gal_region.coord_list[0], gal_region.coord_list[1])
    gal_mask_sum = np.sum(gal_mask)

    # Loop back over the saturation regions, not keeping those that aren't being retained, or which wholly encompass target galaxy
    sat_regions_out = pyregion.ShapeList([])
    for sat_region in sat_regions:

        # Check that saturation region does't wholly encompass target galaxy
        sat_mask = ChrisFuncs.Photom.EllipseMask( np.zeros(image.shape), np.max(sat_region.coord_list[2:4]), np.max(sat_region.coord_list[2:4])/np.min(sat_region.coord_list[2:4]), sat_region.coord_list[4], sat_region.coord_list[0], sat_region.coord_list[1])
        check_mask = gal_mask + sat_mask
        if np.where(check_mask==2)[0].shape[0]<gal_mask_sum:

            # Record regions that pass all criteria
            if int(sat_region.attr[1]['text']) in sat_indices_out:
                sat_regions_out.append(sat_region)

    # Now read in and loop over the star region file, and make new region list using only wanted entries
    star_regions = pyregion.open(star_path)
    star_regions_out = pyregion.ShapeList([])
    for star_region in star_regions:
        if 'text' in star_region.attr[1].keys():
            if int(star_region.attr[1]['text']) in sat_indices_out:
                star_regions_out.append(star_region)
#        elif star_region.name=='point':
#            star_regions_out.remove(star_region)

    # Save updated regions to file
    sat_regions_out.write(sat_path.replace('.reg','_revised.reg'))
    star_regions_out.write(star_path.replace('.reg','_revised.reg'))





# Define function that acquires online catalogues for astromagic processing
def PreCatalogue(source_dict, bands_dict, kwargs_dict):
    if kwargs_dict['verbose']: print '['+source_dict['name']+'] PTS AstroMagic retrieving list of foreground stars in map from online catalogues.'



    # Loop over each band, determining map sizes
    diam_max = 0.0
    for band in bands_dict.keys():
        band_dict = bands_dict[band]

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
        #diam *= 2.0**0.5
        if diam>diam_max:
            diam_max = diam

    # Use IRSA header template service to create dummy header
    dummy_cdelt_arcsec = (diam_max/100.0)*3600.0
    url = 'http://irsa.ipac.caltech.edu/cgi-bin/HdrTemplate/nph-hdr?location='+str(source_dict['ra'])+'%2C+'+str(source_dict['dec'])+'&system=Equatorial&equinox=2000.&width='+str(diam_max)+'&height='+str(diam_max)+'&resolution='+str(dummy_cdelt_arcsec)+'&rotation=0.0'
    sys.stdout = open(os.devnull, "w")
    ChrisFuncs.wgetURL(url, os.path.join(kwargs_dict['temp_dir_path'],'Header_Template.txt'), clobber=True, auto_retry=False)
    sys.stdout = sys.__stdout__
    """
    # Produce dummy header
    dummy_cdelt = diam_max / 11.0
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
    """
    # Create dummy map
    dummy_header = astropy.io.fits.Header.fromfile( os.path.join(kwargs_dict['temp_dir_path'],'Header_Template.txt'), sep='\n', endcard=False, padding=False)
    dummy_map = np.zeros([ int(dummy_header['NAXIS1']), int(dummy_header['NAXIS2']) ])
    dummy_cdelt_arcsec = (diam_max/float(dummy_header['NAXIS2']))*3600.0
    astropy.io.fits.writeto( os.path.join( kwargs_dict['temp_dir_path'],'FITS_Template.fits' ), dummy_map, header=dummy_header)

    # Get AstroMagic catalogue object using dummy fits as reference
    logging.setup_log(level="ERROR")
    importer = ImageImporter()
    importer.run( os.path.join( kwargs_dict['temp_dir_path'],'FITS_Template.fits' ) )
    image = importer.image
    catalog_importer = CatalogImporter()
    catalog_importer.run(image.frames.primary)

    # Return resulting imported catalogue object
    return catalog_importer



