# moments.py
#
# written by: Oliver Cordes 2017-09-26
# changed by: Oliver Cordes 2017-10-11

import numpy as np

class Moments:
    def __init__( self, Amp, x0, y0, Mxx, Mxy, Myy, rho4 ):
        self.moments_amp = Amp
        self.moments_centroid_x = x0
        self.moments_centroid_y = y0
        self.moments_rho4 = rho4
        self.moments_m_xx = Mxx
        self.moments_m_xy = Mxy
        self.moments_m_yy = Myy

        self.shear = Shear( e1=self.e1, e2=self.e2 )

    @property
    def e1( self ):
        return  (self.moments_m_xx-self.moments_m_yy) / (self.moments_m_xx+self.moments_m_yy)

    @property
    def e2( self ):
        return 2.*self.moments_m_xy / (self.moments_m_xx+self.moments_m_yy)

    @property
    def sigma( self ):
        return pow((self.moments_m_xx*self.moments_m_yy-self.moments_m_xy*self.moments_m_xy), 0.25 )

    @property
    def observed_shape( self ):
        return self.shear

    @property
    def moments( self ):
        output_dict = {}
        output_dict['flux'] = self.moments_amp
        output_dict['x'] = self.moments_centroid_x - 0.5 # Compensating for GalSim's default origin
        output_dict['y'] = self.moments_centroid_y - 0.5 # Center of first pixel is at (0.5, 0.5), not (1, 1)
        output_dict['e1'] = self.e1
        output_dict['e2'] = self.e2
        output_dict['sigma'] = self.sigma
        output_dict['rho4'] = self.moments_rho4
        output_dict['m_xx'] = self.moments_m_xx
        output_dict['m_yy'] = self.moments_m_yy
        output_dict['m_xy'] = self.moments_m_xy

        output_dict['g1'] = self.observed_shape.g1
        output_dict['g2'] = self.observed_shape.g2

        return output_dict


class Shear( object ):
    def __init__( self, *args, **kwargs ):

        self._e1 = None   # by default it is not set
        self._e2 = None   # by default it is not set

        # taken from the galsim source code

        if len(kwargs) > 2:
            raise TypeError( "Shear constructor received >2 keyword arguments: %s"%kwargs.keys())

        # Unnamed arg must be a complex shear
        if len(args) == 1:
            self._g = args[0]
            if not isinstance(self._g, complex):
                raise TypeError("Non-keyword argument to Shear must be complex g1 + 1j * g2")

        # Empty constructor means shear == (0,0)
        elif not kwargs:
            self._g = 0j

        elif 'e1' in kwargs or 'e2' in kwargs:
            self._e1 = kwargs.pop( 'e1', 0. )
            self._e2 = kwargs.pop( 'e2', 0. )
            absesq = self._e1**2 + self._e2**2
            if absesq > 1.:
                raise ValueError("Requested distortion exceeds 1: %s"%np.sqrt(absesq))
            self._g = (self._e1 + 1j * self._e2) * self._e2g(absesq)

    @property
    def shear( self ):
        return self._g

    @property
    def g1( self ):
        return self._g.real

    @property
    def g2( self ):
        return self._g.imag

    @property
    def e1( self ):
        if self._e1 is None:
            return self._g.real * self._g2e( abs( self._g )**2 )
        else:
            return self._e1

    @property
    def e2( self ):
        if self._e2 is None:
            return self._g.imag * self._g2e( abs( self._g )**2 )
        else:
            return self._e2

    # helper functions
    def _g2e( self, absgsq):
        return 2. / (1.+absgsq)

    def _e2g(self, absesq):
        if absesq > 1.e-4:
            #return (1. - np.sqrt(1.-absesq)) / absesq
            return 1. / (1. + np.sqrt(1.-absesq))
        else:
            # Avoid numerical issues near e=0 using Taylor expansion
            return 0.5 + absesq*(0.125 + absesq*(0.0625 + absesq*0.0390625))
