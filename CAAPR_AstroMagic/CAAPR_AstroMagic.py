# Import the relevant PTS classes and modules
import sys
import os
sys.path.append( str( os.path.join( os.path.split( os.path.dirname(os.path.abspath(__file__)) )[0], 'CAAPR_AstroMagic', 'PTS') ) )
import signal
import shutil
import pts
import astropy.io.fits
from astropy.units import Unit
from pts.magic.core.image import Image
from pts.magic.misc.imageimporter import ImageImporter
from pts.magic.sources.finder import SourceFinder
from pts.magic.catalog.importer import CatalogImporter
from pts.magic.sources.extractor import SourceExtractor
from pts.core.tools import configuration
from pts.core.tools import logging, time, filesystem
from pts.magic.basics.region import Region





# Define function that wraps the AstroMagic calling procedure
def Magic(pod):
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
        signal.alarm(7200)



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


        """
        # THE STAR, SATURATION AND OTHER-SOURCES REGIONS CAN BE ADJUSTED BY THE USER
        # If the regions are adjusted, they have to be reloaded
        #star_region = Region.from_file(star_region_path)
        #saturation_region = Region.from_file(saturation_region_path)
        #other_region = Region.from_file(other_region_path)
        """


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



