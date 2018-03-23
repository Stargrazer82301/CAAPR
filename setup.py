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
                 install_requires = ['astropy>=2.0',
                                     'aplpy>=1.0',
                                     'colorpy>=0.1.1',
                                     'lmfit>=0.9.9',
                                     'config>=0.3.9',
                                     'pyregion>=2.0',
                                     'photutils>=0.4'],
                 zip_safe = False)
