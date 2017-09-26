# moments.py
#
# written by: Oliver Cordes 2017-09-26
# changed by: Oliver Cordes 2017-09-28


class Moments:
    def __init__( self ):
        self.moments_amp = -1.
        self.moments_centroid_x = -1.
        self.moments_centroid_y = -1.
        self.moments_sigma = -1.
        self.moments_rho4 = -1.
        self.moments_m_xx = -1.
        self.moments_m_yy = -1.
        self.moments_m_xy = -1.

    @property
    def e1( self ):
        return  (self.moments_m_xx-self.moments_m_yy) / (self.moments_m_xx+self.moments_m_yy)

    @property
    def e2( self ):
        return 2.*self.moments_m_xy / (self.moments_m_xx+self.moments_m_yy)

    @property
    def sigma( self ):
        return pow((self.moments_m_xx*self.moments_m_yy-self.moments_m_xy*self.moments_m_xy), 0.24 )

    @property
    def moments( self ):
        output_dict = {}
        output_dict['x'] = self.moments_centroid_x - 0.5 # Compensating for GalSim's default origin
        output_dict['y'] = self.moments_centroid_y - 0.5 # Center of first pixel is at (0.5, 0.5), not (1, 1)
        output_dict['e1'] = self.e1
        output_dict['e2'] = self.e2
        output_dict['sigma'] = self.moments_sigma
        output_dict['rho4'] = self.moments_rho4
        output_dict['m_xx'] = self.moments_m_xx
        output_dict['m_yy'] = self.moments_m_yy
        output_dict['m_xy'] = self.moments_m_xy

        return output_dict
