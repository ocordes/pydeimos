# moments.py
#
# written by: Oliver Cordes 2017-09-26
# changed by: Oliver Cordes 2017-10-02


class Moments:
    def __init__( self, Amp, x0, y0, Mxx, Mxy, Myy, rho4 ):
        self.moments_amp = 2*Amp
        self.moments_centroid_x = x0
        self.moments_centroid_y = y0
        self.moments_rho4 = rho4
        self.moments_m_xx = Mxx
        self.moments_m_xy = Mxy
        self.moments_m_yy = Myy

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

        return output_dict
