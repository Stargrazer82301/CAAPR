# CAAPR (the Comprehensive & Adaptable Aperture Photometry Routine)

CAAPR (the Comprehensive & Adaptable Aperture Photometry Routine) is an aperture photometry pipeline, designed with multiwavelength extragalactic astronomy in mind. It's written in Python 2.7. (Yes, I know that Python 2.7 is *soo* 2016; I'll get around to updating it some day. Meanwhile, if you want to avoid compatibility issues, I recommend using a [virtual environment](https://conda.io/docs/user-guide/tasks/manage-environments.html); it's nice and straightforward.)

CAAPR was created to be able to generate aperture-matched photometry that is truly cross-comparable, even in the face of the broad range of data types that are required for multiwavelength astronomy – producing fluxes ***and uncertainties*** in a consistent manner, despite the enormous variation in the characteristics of observations ranging from ultraviolet to infrared to microwave. CAAPR features a novel method for reliably extrapolating the aperture noise for observations that cover a very limited amount of background.

CAAPR is designed to be an end-to-end photometry pipeline; the user sets it running, and it outputs a couple of CSV tables describing the results, along with some thumbnail images illustrating the apertures. Of course, astronomers tend to distrust "black box" code, so CAAPR is designed with a host of options and switches that the user can adjust, if they want. This makes CAAPR more of a "grey box" pipeline.

CAAPR is described in detail in [C J R Clark et al (2018)](http://adsabs.harvard.edu/abs/2018A%26A...609A..37C), and is a development of the pipeline first used in [C J R Clark et al (2015)](http://adsabs.harvard.edu/abs/2015MNRAS.452..397C). Photometry produced using CAAPR has also been presented in [De Vis et al (2017a)](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1705.02340), [Keenan (2017)](http://adsabs.harvard.edu/abs/2017PhDT........54K), and [Saintonge et al (2018)](http://adsabs.harvard.edu/abs/2018arXiv180907336S).

## Installation

1. [Download](https://github.com/Stargrazer82301/CAAPR/archive/master.zip) the zip archive containing CAAPR from the GitHub repository, and extract it to a temporary location of your choice.
2. Open a terminal and navigate to extracted directory that contains setup.py
3. Install using the terminal command `pip install -e . --user` (Whilst you can try installing using the more common `python setup.py install` route, this sometimes hits compatability issues, whereas the pip local installer seems more robust.)
4. Download and install the ChrisFuncs package, found [here](https://github.com/Stargrazer82301/ChrisFuncs), which CAAPR depends
upon. Like CAAPR it can be downloaded, extracted, and installed very straightforwardly.

## Example Usage

The CAAPR repository has a folder called CAAPR_Example. This contains an example script, and some example data. Once CAAPR has been installed, the example script `CAAPR_Example.py` can be run to provide the user with a basic example of CAAPR in action.

`CAAPR_Example.py` just imports CAAPR,  and then calls the `CAAPR.CAAPR_Main.Run` function which runs the pipeline. The `CAAPR.CAAPR_Main.Run` function has a two required keyword arguments (`sources_table_path` and `bands_table_path`), plus a number of option keyword arguments which allow the user to adjust the pipeline's behavour. The various arguments are as follows:
 - `sources_table_path` [**str, required**] The path (relative or absolute) to a CSV file that provides the details of each target source that CAAPR should process. For column details see subsection below.
 - `bands_table_path` [**str, required**] The path (relative or absolute) to a CSV file that provides the details of each band (ie GALEX FUV, SDSS r, etc) that CAAPR should process. For column details see subsection below.
 - `output_dir_path` [**str**, default = current working directory] The path to a directory where the output tables and thumbnail images produced by CAAPR will be saved. If the directory does not already exist, it will be created.
 - `temp_dir_path` [**str**, default = current working directory] The path to a directory that CAAPR will use for storing temporary working files whilst it is running. If the directory does not already exist, it will be created.
 - `fit_apertures` [**bool**, default = True] If True, CAAPR will fit apertures, and save the apertures to the file named in the `aperture_table_path` kwarg (see below). If False, CAAPR will assume the `aperture_table_path` file is a pre-existing table containing an aperture for each source, and will skip aperture-fitting and progress straight to the photometry phase (assuming the `do_photom` is set to True; see below).
 - `aperture_table_path` [**str/None**, default = within `output_dir_path`] The path to a CSV file that will contain the apertures for each source (see details for fit_apertures kwarg above); if None, a file will be created automatically in the directory given by `output_dir_path`. Alternatively, if you have a set of pre-defined apertures you want CAAPR to do photometry with, this kwarg is used to point to it. For column details see subsection below.
 - `do_photom` [**bool**, default = True] If this is set to False and the `fit_apertures` kwarg is set to True, CAAPR will just fit apertures, and not bother doing any actual photometry.
 - `photom_table_path` [**str/None**, default = within `output_dir_path`] The path to a CSV file where CAAPR will save the measured photometry for each source; if None, a file will be created automatically in the directory given by `output_dir_path`. For column details see subsection below.
 - `expansion_factor` [**float**, default = 1.25] When `fit_apertures = True`, CAAPR's aperture-fitting process will automatically expand the best-fit aperture by this factor, in order to encompass low-SNR emission at the periphery of each target source (see Section 3.4 of C J R Clark et al, 2018).
 - `polysub` [**bool**, default = True] Whether or not CAAPR will attempt to fit and subtract a sky polynomial to remove large-scale background from maps (see Section 3.3 of C J R Clark et al, 2018).
 - `starsub` [**bool**, default = True] Whether or not CAAPR will attempt to remove stars from the input maps (see Section 3.2 of C J R Clark et al, 2018). This can be enabled/disabled for particular bands/sources `bands_table_path` and `sources_table_path` tables.
- `extinction_corr` [**bool**, default = True] Whether or not CAAPR should correct fluxes for foreground Galactic extinction (see Section 3.6 of C J R Clark et al, 2018),
- `parallel` [**bool**, default = True] Whether CAAPR should operate in parallel, processing multiple bands simultaneously for each source. The number of parallel threads employed is given by the `n_proc` kwarg.
- `n_proc` [**int**, default = number of available threads, minus 2] If `parallel = True`, CAAPR will use this many CPU threads when processing each source.
- `thumbnails` [**bool**, default = true] Whether CAAPR should produce thumbnail images to illustrate the photometric apertures used (see Figure 6 of C J R Clark et al, 2018); these images will be placed in the output directory defined by `output_dir_path`.
- `verbose` [**bool**, default = True] Whether CAAPR should print verbose output to the console, to update the user on its progress.
- `messy` [**bool**, default = False] If this is set to True, CAAPR will not delete the contents of the temporary folder (defined by `temp_dir_path`) after each source has been processed. This is useful if the user wishes to inspect the star- and background-subtracted maps. Note that if this is set to True, the thumbnail images produced when `thumbnails ==True` may look a bit ugly, with too much whitespace.

### Sources Table

CAAPR requires a source table, in CSV format (with a one-line header of the column names, described below), in order to operate. When calling CAAPR, the source table is pointed to using the `sources_table_path` kwarg, described above. An example of a sources table can be found in the CAAPR_Example folder of the CAAPR repository. The columns in the table should be as follows:
 - `name` [**str**] The name of the source.
 - `ra` [**float**] The right ascension of the source, in decimal degrees.
 - `dec` [**float**] The declination of the source, in decimal degrees.
 - `aperture_bands_exclude` [**str**] A semicolon-separated list of bands which should not be used for aperture fitting in the case of this source. For example, if know there is a bright SMG right next to a particular source, you might exclude the 350um and 500um bands from aperture fitting; this would be done by entering `SPIRE_350;SPIRE_500` here. This column can be left blank for any source.
 - `photom_bands_exclude` [**str**] A semicolon-separated list of bands which should not be used for photometry in the case of this source. For example, you might have an old optical photographic plate scan image of this source that you want to use for aperture-fitting, but which is no use for photometry. This column can be left blank for any source.
 - `starsub_bands_exclude` [**str**] A semicolon-separated list of bands for which CAAPR should not attempt star subtraction for this source. For example, if the star subtraction seems to perform poorly, and does more harm than good in a particular band, for a particular source. This column can be left blank for any source.
 - `fitting_min_semimaj_arcsec` [**float**] A float describing the minimum semimajor axis of the photometric aperture that should be fitted for a given source. This is in case you have an *a priori* reason to want CAAPR to fit a photometric aperture of at least a certain size for a given source. This column can be left blank for any source.


### Bands Table

CAAPR requires a bands table, in CSV format (with a one-line header of the column names, described below), describing each of the bands that is to be processed. When calling CAAPR, the bands table is pointed to using the `bands_table_path` kwarg, described above. An example of a bands table can be found in the CAAPR_Example folder of the CAAPR repository. The columns in the table should be as follows:
 - `band_name` [**str**] The name of the band. CAAPR assumes that input FITS files are named according to the format `[source_name]_[band_name].fits`. Note that .fits.gz compressed FITS files will also be accepted (and will often be faster).
 - `band_dir` [**str**] The path of the directory where all of the FITS files for this band are stored. Alternatively, if you want to use the same FITS file for every source, this can point towards a specific FITS file (see Section 3.1 of C J R Clark et al, 2018).
 - `use_error_map` [**bool**] Whether error maps are available for this band. At the moment, this functionality is redundant, so it should just be set to False for every entry.
 - `make_cutout` [**bool/float**] If set to False, CAAPR will read in the entire FITS file provided, and use that for the photometry. If set to a float, CAAPR will make a cutout of the FITS file with a diameter equal to this number of arcseconds; if your map if far bigger than is needed for photometry of your source, this makes CAAPR run much faster (see Section 3.1 of C J R Clark et al, 2018).
 - `calib_error` [**str**] The calibration uncertainty of the band in question, in relative units (for example, if a band has a calibration uncertainty of 7%, enter 0.07 here). Whilst CAAPR can determine the measurement uncertainty when performing photometry, there will still be some calibration uncertainty associated with the data. The measurement and calibration uncertainties will be added in quadrature to yield the final uncertainty value it outputs (see Section 3.7 of C J R Clark et al, 2018). If you don't want to do this, just provide a value of 0.0 here.
 - `beam_arcsec` [**float**] The FWHM of this bands' PSF, in arcseconds.
 - `colour_map` [**str**] The colourmap to use for this band when thumbnail image are produced, if the `thumbnails = True` kwarg is set (see above). Will accept the name of any valid MatPlotLib colourmap. For example, just enter `gist_heat` for standard "astronomy orange" coloured images!
 - `annulus_inner` [**float**] The inner semi-major axis of the background annulus used for local sky-subtraction, as a factor of the semi-major axis of the source aperture (see Section 3.6 of C J R Clark et al, 2018). For example, entering 1.25 here would mean that for a source aperture with semi-major axis of 100 arcseconds, the inner semi-major axis of that sources' background annulus would be 125 arcseconds.
 - `annulus_outer` [**float**] The outer semi-major axis of the background annulus used for local sky-subtraction, as a factor of the semi-major axis of the source aperture (see Section 3.6 of C J R Clark et al, 2018). Works the same as for `annulus_inner`, but defines the outer semi-major axis of the background annulus. If you enter 1.601 for this, and 1.25 for `annulus_inner`, then you will have sensibly-sized background annuli that will always have the same area as your source apertures!
 - `subpixel_factor` [**float**] The factor by which the pixels in the input map will be split into smaller sub-pixels (see Section 3.6 of C J R Clark et al, 2018). For example, Herschel-SPIRE 500um maps typically have very large 14 arcsecond pixels. If your target galaxy is only 60 arcseconds across, this is problematic. By entering a value of 10 here, each 14 arcsecond pixel would be split into a 10x10 grid of 1.4 arcsecond pixels, allowing more precise photometry. This will slow things down, however. If you don't want to split a band's pixels into subpixels, just enter a value of 1.0 here.
 - `gain_factor` [**float**] All output photometry will be multiplied by this amount. This is useful if you maps are in different units than what you're actually interested in. If you don't want to alter the output photometry in this way, just enter a value of 1.0 here.
 - `consider_aperture` [**bool**] Whether or not CAAPR should bother considering this band when fitting its master source aperture (see Section 3.4 of C J R Clark et al, 2018). For example, if you have SDSS g-band, r-band, and i-band data, you may as well just set CAAPR to consider g-band and i-band; it's unlikely including r-band will affect the final aperture, and the aperture fitting will run faster.
 - `remove_stars` [**bool**] Whether or not CAAPR should attempt to remove foreground stars in this band (see Section 3.2 of C J R Clark et al, 2018). For example, even if you want to remove stars in most of your bands, there's no point trying for 100um far-infrared data!
 - `downsample_factor` [**float**] The factor by which CAAPR should downsample maps in this band before processing. For example, you may have some Hubble data at a *much* higher resolution than the rest of your maps (like GALEX, SDSS, WISE, Herschel, etc). You have no use for all that extra resolution, and the Hubble maps will take very long to process. But if you enter a value of 5 here for that Hubble band, the data will be downsampled by a factor of 5 (ie, 25 pixels will be combined into 1 during processing), making everything a lot faster. If you don't want to downsample this band, just enter a value of 1.0 here.
 - `beam_correction` [**bool/str**] If False, CAAPR will not attempt to correct fluxes measured in this band for the effect of the instrumental PSF (which will spread some flux out of the source aperture). If set to True, CAAPR will correct the measured fluxes using the beam FWHM provided by the `beam_arcsec` column (see above), and assume that the PSF in this band follows an Airy disc profile. Or, you can enter a file path here, pointing towards a FITS file which contains this band's PSF (see Section 3.6 of C J R Clark et al, 2018).

### Apertures Table

As CAAPR is an aperture photometry pipeline, it unsurprisingly involves apertures. When the `fit_apertures = True` (see above), CAAPR will output a CSV table describing the aperture it determined for each source. Alternatively, if you already have a bunch of apertures you want to do photometry with using CAAPR, you can provide it to the pipeline using the `aperture_table_path` kwarg (see above). Either way, the CSV table should have the following columns (an example table is provided in the CAAPR_Example folder of the repository.  
 - `name` [**str**] Name of the source in question.
 - `semimaj_arcsec` [**float**] Semi-major axis of the source aperture, in arcseconds.
 - `axial_ratio` [**float**] Axial ratio of the source aperture (ie, a/b).
 - `pos_angle` [**float**] Position angle of the source aperture, in degrees.
