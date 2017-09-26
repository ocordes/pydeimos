import pydeimos
import galsim
import numpy as np

import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)





def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    return lambda x,y: height*np.exp(
        -(((float(center_x)-x)/width_x)**2+((float(center_y)-y)/width_y)**2)/2)

# Prepare a coordinate grid
X, Y = np.mgrid[0:100, 0:100]+0.5 # So that the arrays contain the coordinates of the pixel *centers*: first pixel goes from 0.0 to 1.0
# Draw a Gaussian
stamp = gaussian(100.0, 50.0, 50.0, 5.0, 10.0)(X, Y) + np.random.randn(*X.shape)


# Play with the simplemoments code:


stampsize = 100
me = pydeimos.simplemoments.MomentEngine(stampsize)

weights = me.make_weights(func="Gauss", center=None, sigma=20)
#print(weights)

centro = me.measure_first_moments(weights*stamp)

logger.info("Centroid at {}".format(centro))

(qxx, qyy, qxy) = me.measure_second_moments(weights*stamp, center=centro)

(e1, e2, r) = me.compute_ellipse((qxx, qyy, qxy))

logger.info("e1: {}".format(e1))
logger.info("e2: {}".format(e2))
logger.info("r: {}".format(r))



