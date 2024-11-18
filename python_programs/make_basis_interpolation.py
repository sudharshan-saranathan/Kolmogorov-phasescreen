import numpy as np
import sys

from astropy.io import fits
from scipy.interpolate import RectBivariateSpline

def writefits(arr, name):

    hdu = fits.PrimaryHDU(arr)
    hdl = fits.HDUList([hdu])
    hdl.writeto(name)

def main():

    basis    = fits.open('basis.fits')[0].data
    aperture = fits.open('aperture.fits')[0].data
    
    name  = sys.argv[1]
    size  = basis.shape[1]
    size_ = int(sys.argv[2])

    x_old_res = np.linspace(-size/2, size/2, size)
    y_old_res = np.linspace(-size/2, size/2, size)
    
    x_new_res = np.linspace(-size/2, size/2, size_)
    y_new_res = np.linspace(-size/2, size/2, size_)

    abscissae_lr, ordinates_lr = np.meshgrid(x_old_res, y_old_res )
    abscissae_nr, ordinates_nr = np.meshgrid(x_new_res, y_new_res)

    basis_new = np.zeros((basis.shape[0], size_, size_))
    for i in range(basis.shape[0]):
        spline_basis_functions = RectBivariateSpline(x_old_res, y_old_res, basis[i])
        basis_new[i]  = spline_basis_functions.ev(ordinates_nr, abscissae_nr) * aperture

    writefits(basis_new, name)

main()
