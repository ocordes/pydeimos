"""Simple self-coded moments, for playing around

Originally written for MegaLUT's meas.mom.py, but we never developed it further.
"""

import numpy as np


class MomentEngine:
    """Measures moments on square stamps.
    A MomentEngine object allows you to measure moments on many stamps.
    
    The normal procedure:
    - measure first moments using a relatively wide gaussian weight function centered on the stamp center
    - measure central second moments, centered 
    
    """
    def __init__(self, stampsize):
        self.stampsize = stampsize

        self.make_pos()


    def make_pos(self):
        """Prepares x and y position arrays and stamp centers, as we'll keep reusing them.
        """
        (self.xes, self.yes) = np.mgrid[0:self.stampsize, 0:self.stampsize] + 0.5 # We define that the first pixel is centered on (0.5, 0.5)
        (self.xcenter, self.ycenter) = ((self.stampsize/2.0), (self.stampsize/2.0)) # Consistent with this convention.


    def make_weights(self, func="flat", center=None, sigma=None):
        """Builds an array of weights
        """
        if func is "flat":
            weights = np.ones((stampsize, stampsize), dtype=float)

        elif func is  "Gauss":
        
            if center is None:
                center = (self.xcenter, self.ycenter) # Use the stamp center
                
            centerdists    = np.hypot(self.xes - center[0], self.yes - center[1])
            weights = np.exp((-1.0 * centerdists**2) / (2.0 * sigma**2))
            
        else:
            raise RuntimeError("Unknown func")

        assert weights.shape == (self.stampsize, self.stampsize)
        return weights
        

    def measure_first_moments(self, array, sumarray=None):
        """Plain measurement of first moments on a 2D numpy array
        """
        assert array.shape == (self.stampsize, self.stampsize)
        
        if sumarray == None:
            sumarray = np.sum(array)
            
        x = np.sum(self.xes * array) / sumarray
        y = np.sum(self.yes * array) / sumarray

        return (x, y)


    def measure_second_moments(self, array, sumarray=None, center=None):
        """Plain measurement of central second moments on a 2D numpy array.
        """
        assert array.shape == (self.stampsize, self.stampsize)
    
        if sumarray == None:
            sumarray = np.sum(array)
        if center is None:
            center = (self.xcenter, self.ycenter) # Use the stamp center
        
        qxx = np.sum(array * ((self.xes - center[0])**2)) / sumarray
        qyy = np.sum(array * ((self.yes - center[1])**2)) / sumarray
        qxy = np.sum(array * (self.xes - center[0]) * (self.yes - center[1])) / sumarray

        return (qxx, qyy, qxy)


    def compute_ellipse(self, second_moments):
        """Turns a tuple of second moments into (e1, e2, r)
        """
        
        (qxx, qyy, qxy) = second_moments
        
        qdiff = (qxx * qyy) - qxy**2
        qsum = qxx + qyy
        
        if (qdiff < 0.0) or (qsum < 0.0):
            raise ValueError()
        else:
            edenum = qxx + qyy + 2.0*np.sqrt(qdiff)
            e1 = (qxx - qyy) / edenum
            e2 = 2.0 * qxy / edenum
            # and 
            r = np.sqrt(qsum) # inspired by the fact that qxx and qyy are "variances"
            return (e1, e2, r)
            
        
    