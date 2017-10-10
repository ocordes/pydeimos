import logging
logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)

import numpy as np
import astropy.table
import pydeimos


# Preparing a catalog with very simple random galaxy parameters, on a (nstamps_x, nstamps_y) grid
stampsize = 100.0
nstamps_x = 5
nstamps_y = 5

catalog = []
for ix in range(nstamps_x):
    for iy in range(nstamps_y):
        row = {}
        row["ix"] = ix
        row["iy"] = iy
        row["stamp_center_x"] = stampsize/2.0 + ix*stampsize
        row["stamp_center_y"] = stampsize/2.0 + iy*stampsize
        row["tru_gal_x"] = row["stamp_center_x"] + np.random.uniform(-5.0, 5.0)
        row["tru_gal_y"] = row["stamp_center_y"] + np.random.uniform(-5.0, 5.0)
        row["tru_gal_flux"] = np.random.uniform(5000.0, 10000.0)
        row["tru_gal_hlr"] = np.random.uniform(3.0, 6.0)
        row["tru_gal_sersicn"] = np.random.uniform(1.0, 3.0)
        row["tru_gal_g1"] = np.random.uniform(-0.3, 0.3)
        row["tru_gal_g2"] = np.random.uniform(-0.3, 0.3)
        row["tru_psf_sigma"] = np.random.uniform(5.0, 10.0)
        row["tru_psf_g1"] = np.random.uniform(-0.2, 0.2)
        row["tru_psf_g2"] = np.random.uniform(-0.2, 0.2)
        catalog.append(row)
catalog = astropy.table.Table(catalog, masked=True)
catalog.meta.update({"stampsize":stampsize, "nstamps_x":nstamps_x, "nstamps_y":nstamps_y})

print(catalog)

pydeimos.utils.write_pickle(catalog, "run_input_catalog.pkl")
# catalog.write("run_input_catalog.fits") # Would also work here, but it is less general, cannot easily deal with masked arrays.

 

