import pydeimos
import galsim
import numpy as np

import logging
logging.basicConfig(level=logging.DEBUG)





def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    return lambda x,y: height*np.exp(
        -(((float(center_x)-x)/width_x)**2+((float(center_y)-y)/width_y)**2)/2)


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
 

res = galsim.hsm.FindAdaptiveMom(gs_stamp)
   
output_dict = {}
output_dict["flux"] = res.moments_amp
output_dict["x"] = res.moments_centroid.x - 0.5 # Compensating for GalSim's default origin 
output_dict["y"] = res.moments_centroid.y - 0.5 # Center of first pixel is at (0.5, 0.5), not (1, 1)
output_dict["g1"] = res.observed_shape.g1 
output_dict["g2"] = res.observed_shape.g2
output_dict["sigma"] = res.moments_sigma
output_dict["rho4"] = res.moments_rho4   


print(output_dict)