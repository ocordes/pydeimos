
import logging
logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)

import pydeimos
import numpy as np
import astropy.table
import galsim


catalog = pydeimos.utils.read_pickle("run_input_catalog.pkl")
#print(catalog)


stampsize = catalog.meta["stampsize"]
nstamps_x = catalog.meta["nstamps_x"]
nstamps_y = catalog.meta["nstamps_y"]


# Preparing a galsim image to contain the array of stamps
gal_image = galsim.ImageF(stampsize * nstamps_x , stampsize * nstamps_y)
gal_image.scale = 1.0
psf_image = galsim.ImageF(stampsize * nstamps_x , stampsize * nstamps_y)
psf_image.scale = 1.0
obs_image = galsim.ImageF(stampsize * nstamps_x , stampsize * nstamps_y)
obs_image.scale = 1.0


# And loop through the catalog, to draw the stamps
for row in catalog:
                
    # We will draw this galaxy in a postage stamp, but first we need the bounds of this stamp.
    ix = int(row["ix"])
    iy = int(row["iy"])
    bounds = galsim.BoundsI(ix*stampsize+1 , (ix+1)*stampsize, iy*stampsize+1 , (iy+1)*stampsize) # Default Galsim convention, index starts at 1.
    gal_stamp = gal_image[bounds]
    psf_stamp = psf_image[bounds]
    obs_stamp = obs_image[bounds]
    
    # A simple Gaussian profile for the Galaxy
    #gal = galsim.Gaussian(sigma=float(row["tru_gal_sigma"]), flux=float(row["tru_gal_flux"]))
    # A Sersic profile:
    gal = galsim.Sersic(n=row["tru_gal_sersicn"], half_light_radius=row["tru_gal_hlr"], flux=row["tru_gal_flux"])
    
    # We make this profile elliptical
    gal = gal.shear(g1=row["tru_gal_g1"], g2=row["tru_gal_g2"]) # This adds the ellipticity to the galaxy
    # Randomize its position with respect to the stamp center
    gal = gal.shift(row["tru_gal_x"]-row["stamp_center_x"], row["tru_gal_y"]-row["stamp_center_y"])
    # And draw this "true" galaxy on a stamp
    gal.drawImage(gal_stamp)
    
    # We also use a Gaussian for the PSF:
    psf = galsim.Gaussian(flux=1., sigma=row["tru_psf_sigma"])
    psf = psf.shear(g1=row["tru_psf_g1"], g2=row["tru_psf_g2"])
    # And draw this PSF on a stamp:
    psf.drawImage(psf_stamp)
    
    # Convolving the two, and adding noise:
    galconv = galsim.Convolve([gal,psf])
    galconv.drawImage(obs_stamp)
    obs_stamp.addNoise(galsim.GaussianNoise(sigma=1.0))
    

# And saving those images:

print("Writing images...")
gal_image.write("run_gal_image.fits")
psf_image.write("run_psf_image.fits")
obs_image.write("run_obs_image.fits")
print("Done")


 
   