import pydeimos
import galsim
import numpy as np

import logging
logging.basicConfig(level=logging.DEBUG)


from pydeimos.moments import *
from pydeimos.DEIMOS import *


def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    return lambda x,y: height*np.exp(
        -(((float(center_x)-x)/width_x)**2+((float(center_y)-y)/width_y)**2)/2)


# main

# Prepare a coordinate grid
X, Y = np.mgrid[0:120, 0:80]+0.5 # So that the arrays contain the coordinates of the pixel *centers*: first pixel goes from 0.0 to 1.0
assert X[0,0] == 0.5
assert Y.shape == (120, 80)
assert Y[0,79] == 79.5

# Draw a Gaussian
stamp = gaussian(100.0, 65.0, 40.0, 5.0, 10.0)(X, Y) + np.random.randn(*X.shape)
#stamp = gaussian(100.0, 65.0, 40.0, 5.0, 10.0)(X, Y)
# Save it to Fits
pydeimos.utils.write_fits(stamp, "stamp.fits")


print( 'Python FindAdaptiveMom')
hsm = pydeimos.hsm.HSM()
res = hsm.FindAdaptiveMom( stamp.transpose() )

print( res.moments )


psf_moments = Moments( 1.0, 0., 0., 0.0001, 0.0001, 0.1‚, 2.0 )

d = DEIMOS( res )

d.deconvolve( psf_moments )

dres = d.get_moments()

print( dres.moments )


print( 'Done.' )
