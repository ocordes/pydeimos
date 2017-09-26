# moments.py
#
# written by: Oliver Cordes 2017-09-26
# changed by: Oliver Cordes 2017-09-28

from moments import Moments

import numpy  as np
from math import *

class HSM:
    def __init__( self,
                  max_moment_nsig2=25.0,
                  convergence_treshold=1e-6,
                  guess_sig=5. ):
        self.max_moment_nsig2     = max_moment_nsig2
        self.convergence_treshold = convergence_treshold
        self.guess_sig            = guess_sig


    def FindAdaptiveMom( self, object_image ):

        c_x  = object_image.shape[0] / 2.
        c_y  = object_image.shape[1] / 2.
        m_xx = self.guess_sig * self.guess_sig
        m_yy = m_xx
        m_xy = 0.
        self.find_ellipmom_2( object_image,
                              c_x, c_y,
                              m_xx, m_xy, m_yy )
        return Moments()


    def find_ellipmom_2( self,
                         object_image,
                         centroid_x,
                         centroid_y,
                         m_xx,
                         m_xy,
                         m_yy
                         ):

        convergence_factor = 1.0


        Amp = -1000

        while( convergence_factor > self.convergence_treshold ):
            A, Bx, By, Cxx, Cxy, Cyy, rho4 = self.find_ellipmom_1( object_image,
                                                             centroid_x, centroid_y,
                                                             m_xx, m_xy, m_yy )
            break
            pass


    def find_ellipmom_1( self,
                         object_image,
                         x0, y0,
                         Mxx, Mxy, Myy ):


        xmin = 0
        ymin = 0
        xmax = object_image.shape[0]-1
        ymax = object_image.shape[1]-1
        print( 'Entering find_ellipmom_1 with:')
        print( ' x0  : %f' % x0 )
        print( ' y0  : %f' % y0 )
        print( ' Mxx : %f' % Mxx )
        print( ' Mxy : %f' % Mxy )
        print( ' Myy : %f' % Myy )
        print( ' e1  : %f' % ((Mxx-Myy)/(Mxx+Myy)) )
        print( ' e2  : %f' % (2*(Mxy/(Mxx+Myy))) )
        print( ' xmax: %i' % xmax )
        print( ' ymax: %i' % ymax )

        # Compute M^{-1} for use in computing weights
        detM = Mxx * Myy - Mxy*Mxy
        if ( detM <= 0.0) or ( Mxx <= 0. ) or ( Myy<= 0.):
            raise ValueError('non positive definite adaptive moments!' )

        Minv_xx    =  Myy/detM
        TwoMinv_xy = -Mxy/detM * 2.0
        Minv_yy    =  Mxx/detM
        Inv2Minv_xx = 0.5/Minv_xx     # Will be useful later...

        Minv_xx__x_x0__x_x0 = np.array( [ Minv_xx*(x-x0)*(x-x0) for x in range(xmax+1)])

        # Now let's initialize the outputs and then sum
        # over all the pixels
        A = 0.
        Bx = 0.
        By = 0.
        Cxx = 0.
        Cxy = 0.
        Cyy = 0.
        rho4 = 0.

        # rho2 = Minv_xx(x-x0)^2 + 2Minv_xy(x-x0)(y-y0) + Minv_yy(y-y0)^2
        # The minimum/maximum y that have a solution rho2 = max_moment_nsig2 is at:
        #   2*Minv_xx*(x-x0) + 2Minv_xy(y-y0) = 0
        # rho2 = Minv_xx (Minv_xy(y-y0)/Minv_xx)^2 - 2Minv_xy(Minv_xy(y-y0)/Minv_xx)(y-y0)
        #           + Minv_yy(y-y0)^2
        #      = (Minv_xy^2/Minv_xx - 2Minv_xy^2/Minv_xx + Minv_yy) (y-y0)^2
        #      = (Minv_xx Minv_yy - Minv_xy^2)/Minv_xx (y-y0)^2
        #      = (1/detM) / Minv_xx (y-y0)^2
        #      = (1/Myy) (y-y0)^2

        y2 = sqrt(self.max_moment_nsig2 * Myy)  # This still needs the +y0 bit.
        y1 = -y2 + y0
        y2 += y0  #  ok, now it's right.
        iy1 = int( ceil( y1 ) )
        iy2 = int( floor( y2 ) )
        if (iy1 < ymin): iy1 = ymin;
        if (iy2 > ymax): iy2 = ymax;
        print( 'y1,y2  : %f,%f' % (y1,y2) )
        print( 'iy1,iy2: %i,%i' % (iy1,iy2) )

        if ( iy1 > iy2 ):
            raise ValueError( 'bound don\'t make sense')

        for y in range( iy1, iy2+1 ): # for(int y=iy1;y<=iy2;y++)
            y_y0 = y-y0
            TwoMinv_xy__y_y0 = TwoMinv_xy * y_y0
            Minv_yy__y_y0__y_y0 = Minv_yy * y_y0 * y_y0

            # Now for a particular value of y, we want to find the min/max x that satisfy
            # rho2 < max_moment_nsig2.
            #
            # 0 = Minv_xx(x-x0)^2 + 2Minv_xy(x-x0)(y-y0) + Minv_yy(y-y0)^2 - max_moment_nsig2
            # Simple quadratic formula:
            a = Minv_xx
            b = TwoMinv_xy__y_y0
            c = Minv_yy__y_y0__y_y0 - self.max_moment_nsig2
            d = b*b-4.*a*c
            if (d < 0.):
                raise ValueError( 'Failure in finding min/max x for some y!')

            sqrtd = sqrt(d)
            inv2a = Inv2Minv_xx
            x1 = inv2a*(-b - sqrtd) + x0
            x2 = inv2a*(-b + sqrtd) + x0
            ix1 = int( ceil( x1 ) )
            ix2 = int( floor( x2 ) )
            if (ix1 < xmin): ix1 = xmin
            if (ix2 > xmax): ix2 = xmax
            if (ix1 > ix2):
              continue   # rare, but it can happen after the ceil and floor.

            x_x0 = ix1 - x0
            mxxptr = ix1 - xmin

            for x in range( ix1, ix2+1 ):
                # Compute displacement from weight centroid, then
                # get elliptical radius and weight.
                #
                rho2 = Minv_yy__y_y0__y_y0 + TwoMinv_xy__y_y0*x_x0 + Minv_xx__x_x0__x_x0[mxxptr]
                print( 'Using pixel: %i %i with value %f rho2=%f x_x0=%f y_y0=%f' % ( x, y, object_image[x][y], rho2, x_x0, y_y0 ) )

                intensity = exp( -0.5 * rho2 ) * object_image[x][y]

                # Now do the addition
                intensity__x_x0 = intensity * x_x0
                intensity__y_y0 = intensity * y_y0
                A    += intensity
                Bx   += intensity__x_x0
                By   += intensity__y_y0
                Cxx  += intensity__x_x0 * x_x0
                Cxy  += intensity__x_x0 * y_y0
                Cyy  += intensity__y_y0 * y_y0
                rho4+= intensity * rho2 * rho2

                x_x0 += 1

        print( 'exiting find_ellipmom_1:')
        print( ' A    = %f' % A )
        print( ' Bx   = %f' % Bx )
        print( ' By   = %f' % By )
        print( ' Cxx  = %f' % Cxx )
        print( ' Cxy  = %f' % Cxy )
        print( ' Cyy  = %f' % Cyy )
        print( ' rho4 = %f' % rho4 )


        return( A, Bx, By, Cxx, Cxy, Cyy, rho4 )
