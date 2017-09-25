import pydeimos
import galsim
import numpy as np

from math import *

import logging
logging.basicConfig(level=logging.DEBUG)





def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    return lambda x,y: height*np.exp(
        -(((float(center_x)-x)/width_x)**2+((float(center_y)-y)/width_y)**2)/2)


def get_raw_moments( sigma, e1, e2 ):
    sigma_sq = sigma * sigma
    d  = sqrt( 1 - ( e1*e1 + e2*e2 ) )

    m_xx = ( sigma_sq / d ) * ( 1 + e1 )
    m_yy = ( sigma_sq / d ) * ( 1 - e1 )
    m_xy = ( sigma_sq / d ) * e2

    return m_xx, m_yy, m_xy


# Prepare a coordinate grid
X, Y = np.mgrid[0:120, 0:80]+0.5 # So that the arrays contain the coordinates of the pixel *centers*: first pixel goes from 0.0 to 1.0
assert X[0,0] == 0.5
assert Y.shape == (120, 80)
assert Y[0,79] == 79.5

# Draw a Gaussian
stamp = gaussian(100.0, 65.0, 40.0, 5.0, 10.0)(X, Y) + np.random.randn(*X.shape)
# Save it to Fits
pydeimos.utils.write_fits(stamp, "stamp.fits")



# Create a galsim image out of this stamp
gs_stamp = galsim.image.Image(stamp.transpose())
#gs_stamp = galsim.image.Image( stamp.transpose().copy() )

print( 'FindAdaptiveMom')
res = galsim.hsm.FindAdaptiveMom(gs_stamp)

print( dir( res.observed_shape ) )

print( res.observed_shape.getMatrix() )

output_dict = {}
output_dict["flux"] = res.moments_amp
output_dict["x"] = res.moments_centroid.x - 0.5 # Compensating for GalSim's default origin
output_dict["y"] = res.moments_centroid.y - 0.5 # Center of first pixel is at (0.5, 0.5), not (1, 1)
output_dict["g1"] = res.observed_shape.g1
output_dict["g2"] = res.observed_shape.g2
output_dict["e1"] = res.observed_shape.e1
output_dict["e2"] = res.observed_shape.e2
output_dict["sigma"] = res.moments_sigma
output_dict["rho4"] = res.moments_rho4
output_dict["m_xx"] = res.moments_m_xx
output_dict["m_yy"] = res.moments_m_yy
output_dict["m_xy"] = res.moments_m_xy

print( output_dict )

m_xx, m_yy, m_xy = get_raw_moments( res.moments_sigma, res.observed_shape.e1, res.observed_shape.e2 )

output_dict["m_xx_p"] = m_xx
output_dict["m_yy_p"] = m_yy
output_dict["m_xy_p"] = m_xy

print(output_dict)

print( 'Done.' )
