# moments.py
#
# written by: Oliver Cordes 2017-09-26
# changed by: Oliver Cordes 2017-10-02

from moments import Moments

import numpy  as np
from math import *

class HSM:
    def __init__( self,
                  max_moment_nsig2=25.0,
                  convergence_treshold=1e-6,
                  guess_sig=5.,
                  bound_correct_wt=0.25,
                  max_amoment=8000.,
                  max_ashift=15.,
                  max_mom2_iter=400):

        self.max_moment_nsig2     = max_moment_nsig2
        self.convergence_treshold = convergence_treshold
        self.guess_sig            = guess_sig
        self.bound_correct_wt     = bound_correct_wt
        self.max_amoment          = max_amoment
        self.max_ashift           = max_ashift
        self.max_mom2_iter        = max_mom2_iter


    def FindAdaptiveMom( self, object_image ):

        c_x  = ( object_image.shape[1] / 2. ) + 0.5
        c_y  = ( object_image.shape[0] / 2. ) + 0.5
        m_xx = self.guess_sig * self.guess_sig
        m_yy = m_xx
        m_xy = 0.
        Amp, x0, y0, Mxx, Mxy, Myy, rho4 = self.find_ellipmom_2( object_image,
                                                                 c_x, c_y,
                                                                 m_xx, m_xy, m_yy )
        return Moments( Amp, x0, y0, Mxx, Mxy, Myy, rho4 )


    # find_ellipmom_2
    # *** COMPUTES ADAPTIVE ELLIPTICAL MOMENTS OF AN IMAGE ***
    #
    # Finds the best-fit Gaussian:
    #
    # f ~ A / (pi*sqrt det M) * exp( - (r-r0) * M^-1 * (r-r0) )
    #
    # The fourth moment rho4 is also returned.
    # Note that the total image intensity for the Gaussian is 2A.
    #
    # Arguments:
    #   data: ImageView structure containing the image
    # > A: adaptive amplitude
    # > x0: adaptive centroid (x)
    # > y0: adaptive centroid (y)
    # > Mxx: adaptive covariance (xx)
    # > Mxy: adaptive covariance (xy)
    # > Myy: adaptive covariance (yy)
    # Returns:
    #   A: adaptive amplitude
    #   x0: adaptive centroid (x)
    #   y0: adaptive centroid (y)
    #   Mxx: adaptive covariance (xx)
    #   Mxy: adaptive covariance (xy)
    #   Myy: adaptive covariance (yy)
    #   rho4: rho4 moment
    #  convergence_threshold: required accuracy

    def find_ellipmom_2( self,
                         data,
                         x0,
                         y0,
                         Mxx,
                         Mxy,
                         Myy
                         ):

        x00 = x0
        y00 = y0
        num_iter = 0

        convergence_factor = 1.0

        # Set Amp = -1000 as initial value just in case the while() block below is never triggered;
        # in this case we have at least *something* defined to divide by, and for which the output
        # will fairly clearly be junk.
        Amp = -1000

        while( convergence_factor > self.convergence_treshold ):
            # Get moments
            Amp, Bx, By, Cxx, Cxy, Cyy, rho4 = self.find_ellipmom_1( data,
                                                             x0, y0,
                                                             Mxx, Mxy, Myy )

            # Compute configuration of the weight function
            two_psi = atan2( 2* Mxy, Mxx-Myy )
            semi_a2 = 0.5 * ((Mxx+Myy) + (Mxx-Myy)*cos(two_psi)) + Mxy*sin(two_psi)
            semi_b2 = Mxx + Myy - semi_a2

            if (semi_b2 <= 0):
               raise ValueError( 'Error: non positive-definite weight in find_ellipmom_2.' )

            shiftscale = sqrt(semi_b2)
            if (num_iter == 0):
               shiftscale0 = shiftscale

            # Now compute changes to x0, etc.
            dx = 2. * Bx / (Amp * shiftscale)
            dy = 2. * By / (Amp * shiftscale)
            dxx = 4. * (Cxx/Amp - 0.5*Mxx) / semi_b2
            dxy = 4. * (Cxy/Amp - 0.5*Mxy) / semi_b2
            dyy = 4. * (Cyy/Amp - 0.5*Myy) / semi_b2

            #print( 'dx  : %f' % dx )
            #print( 'dy  : %f' % dy )
            #print( 'dxx : %f' % dxx )
            #print( 'dxy : %f' % dxy )
            #print( 'dyy : %f' % dyy )

            if ( dx  >  self.bound_correct_wt ):
                 dx  =  self.bound_correct_wt
            if ( dx  < -self.bound_correct_wt ):
                 dx  = -self.bound_correct_wt
            if ( dy  >  self.bound_correct_wt ):
                 dy  =  self.bound_correct_wt
            if ( dy  < -self.bound_correct_wt ):
                 dy  = -self.bound_correct_wt
            if ( dxx >  self.bound_correct_wt ):
                 dxx =  self.bound_correct_wt
            if ( dxx < -self.bound_correct_wt ):
                 dxx = -self.bound_correct_wt
            if ( dxy >  self.bound_correct_wt ):
                 dxy =  self.bound_correct_wt
            if ( dxy < -self.bound_correct_wt ):
                 dxy = -self.bound_correct_wt
            if ( dyy >  self.bound_correct_wt ):
                 dyy =  self.bound_correct_wt
            if ( dyy < -self.bound_correct_wt ):
                 dyy = -self.bound_correct_wt

            # Convergence tests
            if ( abs( dx ) > abs( dy ) ):
                convergence_factor = abs( dx )
            else:
                convergence_factor = abs( dy )
            convergence_factor *= convergence_factor
            if ( abs( dxx ) > convergence_factor ):
                convergence_factor = abs( dxx )
            if ( abs( dxy ) > convergence_factor ):
                convergence_factor = abs( dxy )
            if ( abs( dyy ) > convergence_factor ):
                convergence_factor = abs( dyy )
            convergence_factor = sqrt( convergence_factor )
            if ( shiftscale<shiftscale0 ):
                convergence_factor *= shiftscale0/shiftscale

            # Now update moments
            x0 += dx * shiftscale
            y0 += dy * shiftscale
            Mxx += dxx * semi_b2
            Mxy += dxy * semi_b2
            Myy += dyy * semi_b2

            # If the moments have gotten too large, or the centroid is out of range,
            #report a failure
            if ( ( abs( Mxx ) > self.max_amoment ) or
                 ( abs( Mxy ) > self.max_amoment ) or
                 ( abs( Myy ) > self.max_amoment ) or
                 ( abs( x0 - x00 ) > self.max_ashift ) or
                 ( abs( y0 - y00 ) > self.max_ashift ) ):
                raise ValueError( 'Error: adaptive moment failed' )

            num_iter += 1
            if ( num_iter > self.max_mom2_iter ):
                raise ValueError( 'Error: too many iterations in adaptive moments' )

        # Re-normalize rho4
        rho4 /= Amp
        return Amp, x0, y0, Mxx, Mxy, Myy, rho4


     # find_ellipmom_1
    # *** FINDS ELLIPTICAL GAUSSIAN MOMENTS OF AN IMAGE ***
    #
    # Returns the parameters:
    # A = int f(x,y) w(x,y)
    # B_i = int (r_i-r0_i) f(r) w(r)
    # C_ij = int (r_i-r0_i) (r_j-r0_j) f(r) w(r)
    # rho4 = int rho^4 f(r) w(r)
    #
    # where w(r) = exp(-rho^2/2), rho^2 = (x-x0) * M^{-1} * (y-y0),
    # M = adaptive covariance matrix, and note that the weight may be set to zero for rho^2 >
    # hsmparams->max_moment_nsig2 if that parameter is defined.
    #
    # Arguments:
    #   data: the input image (ImageView format)
    #   x0: weight centroid (x coordinate)
    #   y0: weight centroid (y coordinate)
    #   Mxx: xx element of adaptive covariance
    #   Mxy: xy element of adaptive covariance
    #  Myy: yy element of adaptive covariance
    # Returns:
    # > A: amplitude
    # > Bx: weighted centroid displacement (x)
    # > By: weighted centroid displacement (y)
    # > Cxx: weighted covariance (xx)
    # > Cxy: weighted covariance (xy)
    # > Cyy: weighted covariance (yy)
    # > rho4w: weighted radial fourth moment

    def find_ellipmom_1( self,
                         object_image,
                         x0, y0,
                         Mxx, Mxy, Myy ):


        xmin = 1
        ymin = 1
        xmax = object_image.shape[1]
        ymax = object_image.shape[0]
        #print( 'Entering find_ellipmom_1 with:')
        #print( ' x0  : %f' % x0 )
        #print( ' y0  : %f' % y0 )
        #print( ' Mxx : %f' % Mxx )
        #print( ' Mxy : %f' % Mxy )
        #print( ' Myy : %f' % Myy )
        #print( ' e1  : %f' % ((Mxx-Myy)/(Mxx+Myy)) )
        #print( ' e2  : %f' % (2*(Mxy/(Mxx+Myy))) )
        #print( ' xmax: %i' % xmax )
        #print( ' ymax: %i' % ymax )

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
        #print( 'y1,y2  : %f,%f' % (y1,y2) )
        #print( 'iy1,iy2: %i,%i' % (iy1,iy2) )

        if ( iy1 > iy2 ):
            raise ValueError( 'bound don\'t make sense')

        for y in range( iy1, iy2+1 ): # for(int y=iy1;y<=iy2;y++)
            y_y0 = y-y0
            TwoMinv_xy__y_y0 = TwoMinv_xy * y_y0
            Minv_yy__y_y0__y_y0 = Minv_yy * y_y0 * y_y0

            #print( 'TwoMinv_xy__y_y0   : %f' % TwoMinv_xy__y_y0 )
            #print( 'Minv_vy__y_y0__y_y0: %f' % Minv_yy__y_y0__y_y0 )

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
                mxxptr += 1
                rho2 = Minv_yy__y_y0__y_y0 + TwoMinv_xy__y_y0*x_x0 + Minv_xx__x_x0__x_x0[mxxptr]
                #print( 'Using pixel: %i %i with value %f rho2=%f x_x0=%f y_y0=%f' % ( x, y, object_image[y-1][x-1], rho2, x_x0, y_y0 ) )
                intensity = exp( -0.5 * rho2 ) * object_image[y-1][x-1]

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

        #print( 'exiting find_ellipmom_1:')
        #print( ' A    = %f' % A )
        #print( ' Bx   = %f' % Bx )
        #print( ' By   = %f' % By )
        #print( ' Cxx  = %f' % Cxx )
        #print( ' Cxy  = %f' % Cxy )
        #print( ' Cyy  = %f' % Cyy )
        #print( ' rho4 = %f' % rho4 )


        return( A, Bx, By, Cxx, Cxy, Cyy, rho4 )
