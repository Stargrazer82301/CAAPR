'''
Hard link the processed CAAPR maps to a more standardized directory structure,
where each galaxy has its own directory.
'''

import argparse
import os

parser = argparse.ArgumentParser(description='Link processed CAAPR maps to a per-galaxy directory structure.')
parser.add_argument('--runname', type=str, default='dustpedia_1')
# Place in data_basedir/galname/target_dir/bandname.fits
# E.g. if data_basedir is ./dustpedia/ (but probably should run from within dustpedia):
# ./dustpedia/NGC1365/native/GALEX_NUV.fits
parser.add_argument('--data_basedir', type=str, default=os.path.abspath('./'))
parser.add_argument('--target_dir', type=str, default='native')
args = parser.parse_args()

caapr_mapsdir = os.path.join(os.path.expanduser('~/CAAPR_run/'), args.runname, 'CAAPR_Output/Processed_Maps/')
data_basedir = args.data_basedir
target_dir = args.target_dir
image_names = list(filter(lambda x: x.endswith('.fits'), os.listdir(caapr_mapsdir)))
n_images = len(image_names)

if n_images == 0:
    raise ValueError("No images found in " + caapr_mapsdir + "!")
print('Linking {} galaxies from {} to {}...'.format(n_images, caapr_mapsdir, data_basedir))

for imgname in image_names:
    # 1 split: ESO123-455_GALEX_NUV.fits -> ESO123-455, GALEX_NUV.fits
    galname, bandname = imgname.split('_', 1)
    srcpath = os.path.join(caapr_mapsdir, imgname)
    target_dirpath = os.path.join(data_basedir, galname, target_dir)
    if not os.path.exists(target_dirpath):
        os.makedirs(target_dirpath)
    targetpath = os.path.join(target_dirpath, bandname)
    os.link(srcpath, targetpath)