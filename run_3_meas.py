
import logging
logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)

import pydeimos
import numpy as np
import astropy.table
import galsim



catalog = pydeimos.utils.read_pickle("run_input_catalog.pkl")
#print(catalog)

# We add masked empty columns for the measurements:
for colname in ["hsm_gal_flux", "hsm_gal_x", "hsm_gal_y", "hsm_gal_g1", "hsm_gal_g2",
                "hsm_gal_sigma", "hsm_gal_rho4",
                "hsm_ksb_g1", "hsm_ksb_g2"]:
    catalog.add_column(astropy.table.MaskedColumn(name=colname, mask=True, dtype=float, length=len(catalog)))

stampsize = catalog.meta["stampsize"]
nstamps_x = catalog.meta["nstamps_x"]
nstamps_y = catalog.meta["nstamps_y"]

psf_image = galsim.fits.read("run_psf_image.fits")
obs_image = galsim.fits.read("run_obs_image.fits")


# And loop through the catalog, to measure 
for row in catalog:
                
    # We will draw this galaxy in a postage stamp, but first we need the bounds of this stamp.
    ix = int(row["ix"])
    iy = int(row["iy"])
    bounds = galsim.BoundsI(ix*stampsize+1 , (ix+1)*stampsize, iy*stampsize+1 , (iy+1)*stampsize) # Default Galsim convention, index starts at 1.
    psf_stamp = psf_image[bounds]
    obs_stamp = obs_image[bounds]
    
    
    # We measure this stamp with HSM
    try:
        res = galsim.hsm.FindAdaptiveMom(obs_stamp, guess_sig=15.0)
    except:
        print("Failed on stamp ({},{})".format(row["ix"], row["iy"]))
        continue
    row["hsm_gal_flux"] = res.moments_amp
    row["hsm_gal_g1"] = res.observed_shape.g1
    row["hsm_gal_g2"] = res.observed_shape.g2
    row["hsm_gal_x"] = res.moments_centroid.x -0.5
    row["hsm_gal_y"] = res.moments_centroid.y -0.5
    row["hsm_gal_sigma"] = res.moments_sigma
    row["hsm_gal_rho4"] = res.moments_rho4
    
    # Let's run a KSB shear estimator on this:
    res = galsim.hsm.EstimateShear(obs_stamp, psf_stamp, shear_est="KSB")
    row["hsm_ksb_g1"] = res.corrected_g1
    row["hsm_ksb_g2"] = res.corrected_g2


# And save the catalog:
#print(catalog)
pydeimos.utils.write_pickle(catalog, "run_meas_catalog.pkl")





 
   