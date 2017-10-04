import pydeimos
import galsim
import numpy as np

import logging
logging.basicConfig(level=logging.DEBUG)


import hsm



def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    return lambda x,y: height*np.exp(
        -(((float(center_x)-x)/width_x)**2+((float(center_y)-y)/width_y)**2)/2)



# check galsim/python hsm outputs
def check_galsim_pyhsm( galsim_dict, pyhsm_dict ):
    print( 'check galsim/pyhsm (diffs):' )
    for i in ['flux', 'sigma', 'x', 'y', 'm_xx', 'm_xy', 'm_yy', 'rho4', 'e1', 'e2' ]:
        print( ' %-10s: %20.15f' % ( i, ( galsim_dict[i] - pyhsm_dict[i] ) ) )


def generate_stamp( size_x=120, size_y=80,
                    x=65., y=40., height=100., width_x=5., width_y=10., angle=0.,
                    have_noise=True ):
    # Prepare a coordinate grid
    X, Y = np.mgrid[0:size_x, 0:size_y]+0.5 # So that the arrays contain the coordinates of the pixel *centers*: first pixel goes from 0.0 to 1.0
    assert X[0,0] == 0.5
    assert Y.shape == (size_x, size_y)
    assert Y[0,size_y-1] == size_y-0.5

    stamp = gaussian(100.0, 65.0, 40.0, 5.0, 10.0)(X, Y)
    if have_noise:
        stamp += np.random.randn(*X.shape)

    return stamp

# main

# generate test variables
vtgen = pydeimos.vtgenerator.VTGenerator()
vtgen.add_value( 'x', 'FLOAT', 60., 70., 1. )
vtgen.add_value( 'y', 'FLOAT', 35., 45., 1. )

for i in vtgen.product_dict():
    stamp = generate_stamp( **i )

stamp = generate_stamp()

# Save it to Fits
pydeimos.utils.write_fits(stamp, "stamp.fits")


# Create a galsim image out of this stamp
gs_stamp = galsim.image.Image(stamp.transpose())


print( 'FindAdaptiveMom')
res = galsim.hsm.FindAdaptiveMom(gs_stamp)

output_dict = {}
output_dict["flux"] = res.moments_amp
output_dict["x"] = res.moments_centroid.x - 0.5 # Compensating for GalSim's default origin
output_dict["y"] = res.moments_centroid.y - 0.5 # Center of first pixel is at (0.5, 0.5), not (1, 1)
#output_dict["g1"] = res.observed_shape.g1
#output_dict["g2"] = res.observed_shape.g2
output_dict["e1"] = res.observed_shape.e1;
output_dict["e2"] = res.observed_shape.e2;
output_dict["sigma"] = res.moments_sigma
output_dict["rho4"] = res.moments_rho4
output_dict["m_xx"] = res.moments_m_xx
output_dict["m_yy"] = res.moments_m_yy
output_dict["m_xy"] = res.moments_m_xy

print(output_dict)

print( 'Python FindAdaptiveMom')
hsm = hsm.hsm.HSM()
res = hsm.FindAdaptiveMom( stamp.transpose() )

output_dict2 = res.moments

print( output_dict2 )


check_galsim_pyhsm( output_dict, output_dict2 )

print( 'Done.' )
