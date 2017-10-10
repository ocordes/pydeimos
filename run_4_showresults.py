

import galsim
import pydeimos
import pydeimos.sky_image_plot as sip


catalog = pydeimos.utils.read_pickle("run_meas_catalog.pkl")
print(catalog)


obs_image = galsim.fits.read("run_obs_image.fits")


# And we visualize the grid and the measurements:
sf = sip.SimpleFigure(obs_image.array, z1=-3.0, z2=10.0, scale=2) # There is no "tranpose" here, as we've defined the stamps with our own convention, and consistently drawn and measured them with galsim.

# The true PSF shape, in blue:
sf.draw_g_ellipses(catalog, x="tru_gal_x", y="tru_gal_y", g1="tru_psf_g1", g2="tru_psf_g2", sigma="tru_psf_sigma", edgecolor="blue", linewidth=1.0)
# The observed galaxy shape, in red:
sf.draw_g_ellipses(catalog, x="hsm_gal_x", y="hsm_gal_y", g1="hsm_gal_g1", g2="hsm_gal_g2", sigma="hsm_gal_sigma", edgecolor="red", linewidth=1.0)


# The true galaxy shape, in thick yellow:
sf.draw_g_ellipses(catalog, x="tru_gal_x", y="tru_gal_y", g1="tru_gal_g1", g2="tru_gal_g2", sigma="tru_gal_hlr", edgecolor="yellow", linewidth=3.0)
# The KSB PSF-corrected shear, in purple:
sf.draw_g_ellipses(catalog, x="tru_gal_x", y="tru_gal_y", g1="hsm_ksb_g1", g2="hsm_ksb_g2", sigma="tru_gal_hlr", edgecolor="purple", linewidth=1.0)


sf.annotate(catalog, x="hsm_gal_x", y="hsm_gal_y", text="g1 = {row[hsm_gal_g1]:.2f}", color="white", xytext=(20, -10))
sf.annotate(catalog, x="hsm_gal_x", y="hsm_gal_y", text="g2 = {row[hsm_gal_g2]:.2f}", color="white", xytext=(20, -25))


sf.show()
