# DEIMOS.py
#
# written by: Oliver Cordes 2017-10-09
# changed by: Oliver Cordes 2017-10-09
#

import moments

class DEIMOS( object ):
    def __init__( self, moments ):
        self.moments = moments

    def get_moments( self ):
        return self.moments

    def deconvolve( self, mpsf ):
        # use explicit relations for up to 2nd moments
        # g(0,0) /= p(0,0);
        self.moments.moments_amp /= mpsf.moments_amp

        #if (Nmin >= 1) {
        #g(0,1) -= g(0,0)*p(0,1);
        self.moments.moments_centroid_y -= self.moments.moments_amp * mpsf.moments_centroid_y
        #g(0,1) /= p(0,0);
        self.moments.moments_centroid_y /= mpsf.moments_amp
        #g(1,0) -= g(0,0)*p(1,0);
        self.moments.moments_centroid_x -= self.moments.moments_amp * mpsf.moments_centroid_x
        #g(1,0) /= p(0,0);
        self.moments.moments_centroid_x /= mpsf.moments_amp

        #if (Nmin >= 2) {
        #g(0,2) -= g(0,0)*p(0,2) + 2*g(0,1)*p(0,1);
        self.moments.moments_m_yy -= ( self.moments.moments_amp * mpsf.moments_m_yy +
                                     2.*self.moments.moments_centroid_y * mpsf.moments_centroid_y )
        #g(0,2) /= p(0,0);
        self.moments.moments_m_yy /= mpsf.moments_amp
        #g(1,1) -= g(0,0)*p(1,1) + g(0,1)*p(1,0) + g(1,0)*p(0,1);
        self.moments.moments_m_xy -= ( self.moments.moments_amp * mpsf.moments_m_xy +
                                    self.moments.moments_centroid_y * mpsf.moments_centroid_y +
                                    self.moments.moments_centroid_x * mpsf.moments_centroid_x )
        #g(1,1) /= p(0,0);
        self.moments.moments_m_xy /= mpsf.moments_amp
        #g(2,0) -= g(0,0)*p(2,0) + 2*g(1,0)*p(1,0);
        self.moments.moments_m_xx -= ( self.moments.moments_amp * mpsf.moments_m_xx +
                                    2.*self.moments.moments_centroid_x * mpsf.moments_centroid_x )
        #g(2,0) /= p(0,0);
        self.moments.moments_m_xx /= mpsf.moments_amp


    def convolve( self, psf_moments ):
        pass
