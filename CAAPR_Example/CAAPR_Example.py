# Import smorgasbord
import multiprocessing as mp
import CAAPR

# Call function which launches CAAPR pipeline
CAAPR.CAAPR_Main.Run(bands_table_path = '../CAAPR_Example/CAAPR_Band_Table_Test.csv', # This points towards a CSV file with a table describing each band that is to be processed
                     sources_table_path = '../CAAPR_Example/CAAPR_Source_Table_Test.csv', # This points towards a CSV file with a table describing each source that is to be processed
                     output_dir_path = '../CAAPR_Example/CAAPR_Output', # This points towards a directory where the output tables and thumbnail images will be saved
                     temp_dir_path = '../CAAPR_Example/CAAPR_Temp', # This points towards a directory that will be used as temporary storage area for the pipeline
                     fit_apertures = True, # If True, CAAPR will fit apertures, and save the apertures to the file named in the aperture_table_path kwarg (below). If False, CAAPR will assume the aperture_table_path file is a pre-existing table containing an aperture for each source, and will skip aperture-fitting and progress straight to the photometry phase
                     aperture_table_path = None, # This points towards a CSV file that will contain the apertures for each source (see details for fit_apertures kwarg above); if None, a file will be
                     photom_table_path = None, # This points towards a file  where CAAPR will save the measured photometry for each soruce
                     expansion_factor = 1.25, # The expansion factor to be applied to aperture during fitting (see Section 3.4 of CJR Clark et al, 2018)
                     polysub = True, # Whether or not CAAPR should attempt to fit and subtract a polynomial to remove background from maps (see Section 3.3 of CJR Clark et al, 2018)
                     starsub = True, # Whether or not CAAPR will attempt to remove stars from bands where it is requested by the bands_table_path table (see Section 3.2 of CJR Clark et al, 2018)
                     do_photom = True, # If this is False and fit_apertures is True, CAAPR will just fit apertures, and not bother doing any actual photometry
                     extinction_corr = True, # Whether or not CAAPR should correct fluxes for foregound Galactic extinction (see Section 3.6 of CJR Clark et al, 2018)
                     parallel = False, # Whether CAAPR should operate in paralle, processing multiple bands simultaneously
                     n_proc = mp.cpu_count()-2, # If parallel is set to True, CAAPR will use this many CPU threads
                     thumbnails = True, # If True, CAAPR will produce thumbnail images to illustrate the results of the photometry (see Figure 6 of CJR Clark et al, 2018)
                     verbose = True, # If True, CAAPR will print a lot of text to the console to update the user on its progress; if False, it won't
                     )

# Jubilate
print('All done!')