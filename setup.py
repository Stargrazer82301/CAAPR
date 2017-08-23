# Import smorgasbord
import setuptools

# Configure setup
setuptools.setup(name = 'caapr',
                 version = '1.0',
                 description = 'CAAPR: the Comprehensive & Adaptable Aperture Photometry Routine. A highly-customisable astronomical aperture photometry pipeline.',
                 url = 'https://github.com/Stargrazer82301/CAAPR',
                 author = 'Chris Clark (github.com/Stargrazer82301)',
                 author_email = 'cjrc88@gmail.com',
                 license = 'MIT',
                 classifiers = ['Programming Language :: Python :: 2.7',
                                'License :: OSI Approved :: MIT License'],
                 packages = setuptools.find_packages(),
                 #packages = ['CAAPR_AstroMagic.pts'].append(setuptools.find_packages()),
                 #package_dir = {'pts':os.path.join('CAAPR_AstroMagic','PTS','pts')},
                 zip_safe = False)
