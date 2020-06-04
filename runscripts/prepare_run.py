'''
Creates directory structure for run, and fetches galaxies.

Arguments:
- galaxies : a single column csv file, containing the dustpedia galaxies that will be processed. Defaults to './galaxy_list.csv'. The directory in which this file is placed will act as the base directory for the run.
- caaprdir : the directory where CAAPR can be found. Default: '~/CAAPR/'.
'''

import argparse
import os
import time
import shutil
import urllib
from functools import partial
import astropy.io.fits
import numpy as np

def verify_exists(path):
    if not os.path.exists(path):
        raise IOError("Path " + path + " does not exist")

def get_basedir(galaxy_filename):
    return os.path.dirname(os.path.abspath(galaxy_filename))

def create_dirstructure(basedir):
    """Creates the necessary directories but leaves them empty."""
    
    dirs = ['tables', 'Cutouts']
    for dirname in dirs:
        if not os.path.exists(dirname):
            os.mkdir(dirname)

def copy_configfiles(basedir, caaprdir):
    """Copy config files from caaprdir to basedir."""

    def copyfile(filename, caapr_subdir='DustPedia', base_subdir=''):
        frompath = os.path.join(caaprdir, caapr_subdir, filename)
        shutil.copy(frompath, os.path.join(basedir, base_subdir))

    copyfile('runconfig.json')
    copyfile('CAAPR_Aperture_Table.csv', base_subdir='tables')
    copyfile('CAAPR_Band_Table.csv', base_subdir='tables')
    copyfile('CAAPR_Source_Table.csv', base_subdir='tables')
    copyfile('run.py', caapr_subdir='runscripts')

def galaxies_from_file(galaxy_filename):
    return list(np.loadtxt(galaxy_filename, dtype=str))

def subset_configfiles(basedir, galaxies):
    """Discard the galaxies that are not needed from the config tables."""

    def subset_configfile(name):
        filename = os.path.join(basedir, 'tables', 'CAAPR_{}_Table.csv'.format(name))
        # Read lines to list
        with open(filename, 'r') as cfgfile:
            lines = cfgfile.readlines()
        # Filter
        def in_sample(line):
            for gal in galaxies:
                if gal in line:
                    return True
            return False

        lines = [lines[0]] + filter(in_sample, lines[1:])
        # Write to same config file
        with open(filename, 'w') as cfgfile:
            cfgfile.writelines(lines) 
    
    subset_configfile('Aperture')
    subset_configfile('Source')

def log_verbose(msg, verbose):
    """Log message, only if verbose is True."""

    if verbose:
        print(msg)

def fetch_galaxies(basedir, galaxies, verbose=False):
    """Check if the data for each galaxy is present. If not: download from dustpedia archive."""

    log = partial(log_verbose, verbose=verbose)

    # Get bands from band table
    bands_filename = os.path.join(basedir, 'tables/CAAPR_Band_Table.csv')
    bandsdata = np.loadtxt(bands_filename, dtype=str, delimiter=',', usecols=(0, 1))[1:, :]
    bands, banddirs = bandsdata[:, 0], bandsdata[:, 1]
    band_dirmap = {bands[i]: banddirs[i] for i in range(len(bands))}

    for banddir in banddirs:
        if not os.path.exists(banddir):
            os.makedirs(banddir)

    # Download bands
    for galaxy in galaxies:
        for band in bands:
            banddir = band_dirmap[band]
            target_filename = os.path.join(banddir, '{}_{}.fits'.format(galaxy, band))
            missing_filename = os.path.join(banddir, '{}_{}_missing'.format(galaxy, band))
            tele_name = os.path.dirname(banddir).split('/')[-1]
            source_url = ('http://dustpedia.astro.noa.gr/Data/GetImage?imageName={}_{}.fits&instrument={}'
                                .format(galaxy, band, tele_name))
            should_download = True
            if os.path.exists(target_filename):
                # file exists: check if it is readable
                try:
                    astropy.io.fits.getdata(target_filename, header=True)
                    should_download = False
                    log('{} - {} already found.'.format(galaxy, band))
                except (IOError, TypeError):
                    log('{} - {} was found but can not be read. Redownloading...'.format(galaxy, band))
                    os.remove(target_filename)
            elif os.path.exists(missing_filename):  # Band previously found to be missing
                log('{} - {} previously found to not be in database.'.format(galaxy, band))
                should_download = False
            if should_download:
                log('{} - downloading {}...'.format(galaxy, band))
                code, content = download_retry(source_url, verbose=verbose)
                if code == 500:
                    log('{} - {} not in database.'.format(galaxy, band))
                    open(missing_filename, 'w').close()  # Create empty file to signal missing band
                elif code != 200:
                    raise IOError("Can not download {} of {}, but it should be there (HTTP code {}).".format(band, galaxy, code))
                else:
                    with open(target_filename, 'wb') as target_file:
                        target_file.write(content)

def download_retry(source_url, max_retries=100, verbose=False):
    log = partial(log_verbose, verbose=verbose)

    i = 1
    while i <= max_retries:
        try:
            response = urllib.urlopen(source_url)
            code = response.getcode()
            if code == 200:
                content = response.read()
            else:
                content = None
            return code, content
        except IOError as e:
            log('Connection failed ({}): '.format(i) + str(e))

        i += 1
        time.sleep(1)
    raise IOError("Can not establish database connection after {} tries.".format(max_retries))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare for CAAPR run: create directory structure and fetch galaxies.')
    parser.add_argument('--galaxies', type=str, default='./galaxy_list.csv')
    parser.add_argument('--caaprdir', type=str, default=None)
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()
    
    galaxy_filename = args.galaxies
    caaprdir = args.caaprdir
    if caaprdir is None:
        # If the script is run from caaprdir/runscripts/, then we just need to go up one directory
        caaprdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        if args.verbose:
            print('Using '+caaprdir+' as CAAPR directory.')
    verify_exists(galaxy_filename)
    verify_exists(caaprdir)
    basedir = get_basedir(galaxy_filename)
    os.chdir(basedir)
    galaxies = np.loadtxt(galaxy_filename, dtype=str)
    if args.verbose:
        print('Preparing {} galaxies.'.format(len(galaxies)))

    create_dirstructure(basedir)
    copy_configfiles(basedir, caaprdir)
    subset_configfiles(basedir, galaxies)
    fetch_galaxies(basedir, galaxies, args.verbose)