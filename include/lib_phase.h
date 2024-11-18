#ifndef _SIMULATOR_
#define _SIMULATOR_

#include <chrono>
#include <random>
#include <cstring>
#include <complex>

#include "fftw3.h"
#include "config.h"
#include "lib_array.h"

#define MAX(a, b) a > b ? a : b
#define MIN(a, b) a < b ? a : b

#define _USE_MATH_DEFINES

typedef std::chrono::high_resolution_clock Time;

template <class type> 
void make_aperture_function(Array<type>& aperture, type aperture_radius){

/*
 * Vector declaration.
 * --------------------------------------------
 * Name             Type            Description
 * --------------------------------------------
 * aperture_dims    sizt_vector     Dimensions of the aperture
 */

    sizt_vector aperture_dims = aperture.get_dims();

/*
 * Variable declaration.
 * ----------------------------------------
 * Name                 Type    Description
 * ----------------------------------------
 * aperture_center_x    type    Center of aperture, abscissa.
 * aperture_center_y    type    Center of aperture, ordinate.
 * distance_center      type    Distance to the center of aperture.
 */ 

    type aperture_center_x = type(aperture_dims[0])/2.0;
    type aperture_center_y = type(aperture_dims[1])/2.0;
    type distance_center   = 0.0;

/* -----------------------
 * Make aperture function.
 * -----------------------
 */

    for(sizt xs = 0; xs < aperture_dims[0]; xs++){
        for(sizt ys = 0; ys < aperture_dims[1]; ys++){
            distance_center = sqrt(pow(xs - aperture_center_x, 2) + pow(ys - aperture_center_y, 2));
            aperture(xs, ys) = distance_center <= aperture_radius ? static_cast<type>(1) : static_cast<type>(0);
        }
    }
}

void make_phase_screen_fourier_shifted(Array<cmpx>& fourier, precision fried, precision sim_size){

/*
 * Variable declaration
 * ----------------------------------------
 * Name             Type        Description
 * ----------------------------------------
 * time_stamp       Time      System clock.  
 */

    Time::time_point time_start = Time::now();

/*
 * Vector declaration
 * --------------------------------------------
 * Name             Type            Description
 * --------------------------------------------
 * fourier_dims     sizt_vector     Dimensions of the fourier array.
 */
    
    sizt_vector fourier_dims = fourier.get_dims();

/*
 * Variable declaration
 * --------------------------------
 * Name     Type        Description
 * --------------------------------
 * xc       precision   Center of the fourier array.
 * yc       precision   Center of the fourier array.
 * amp      precision   Amplitude of the fourier array.
 * frq      precision   Spatial frequency.
 */

    precision xc  = fourier_dims[0]/2;
    precision yc  = fourier_dims[1]/2;
    precision amp = 1.0;
    precision frq = 1.0;

/*
 * Variable declaration
 * ------------------------------------------------
 * Name             Type                Description
 * ------------------------------------------------
 * duration         Time::duration     Time duration until execution of this line. 
 * mt19937_seed     sizt                Seed for mersenne twister random number generator. 
 */
    
    Time::duration duration = Time::now() - time_start;
    sizt mt19937_seed = duration.count();

/* ----------------------------------------------------
 * Initialize mersenne twister random number generator.
 * ----------------------------------------------------
 */

    std::mt19937 generator(mt19937_seed);
    std::normal_distribution<precision> distribution(0.0, 1.0);

/* --------------------------------------
 * Loop to construct Kolmogorov spectrum.
 * --------------------------------------
 */

    for(sizt xpix = 0; xpix < fourier_dims[0]; xpix++){
        for(sizt ypix = 0; ypix < fourier_dims[1]; ypix++){

        /* --------------
         * Set frequency.
         * --------------
         */

            frq = sqrt(pow(precision(xpix) - xc, 2) + pow(precision(ypix) - yc, 2));
           
        /* --------------
         * Set amplitude.
         * --------------
         */
        
            amp = frq == 0.0 ? 0.0 : sqrt(0.023) * pow(2.0*sim_size/fried, 5./6.) / pow(frq, 11./6.);

        /* 
         * Variable declaration
         * --------------------------------
         * Name     Type        Description
         * --------------------------------
         * cosphi   precision   Real part of the fourier phase.
         * sinphi   precision   Imag part of the fourier phase.
         * phase    cmpx        Fourier transformed phase.
         */
            
            precision cos_phi = distribution(generator);
            precision sin_phi = distribution(generator);
            
            cmpx phase(amp*cos_phi, amp*sin_phi);

        /* ------------------------------------------------------------------------
         * Set the value of the fourier transform at (xpix, ypix), shifted to edge.
         * ------------------------------------------------------------------------
         */

            fourier(sizt(xpix + fourier_dims[0] - xc) % fourier_dims[0], sizt(ypix + fourier_dims[1] - yc) % fourier_dims[1]) = phase; 
        }
    }
}

void make_residual_phase_screen(Array<precision>& phase, Array<precision>& basis, Array<precision>& weights, sizt_vector norm){

/*
 * Vector declaration
 * --------------------------------------------
 * Name             Type            Description
 * --------------------------------------------
 * dims_basis       sizt_vector     Dimensions of the array storing the basis functions.
 * dims_weights     sizt_vector     Dimensions of the array storing the basis weights.
 * dims_mode        sizt_vector     Dimensions of the array storing a single basis mode.
 * dims_mode_cfs    sizt_vector     Dimensions of the array storing the mode amplitudes.
 */

    sizt_vector dims_basis = basis.get_dims();
    sizt_vector dims_weights = weights.get_dims();
    sizt_vector dims_mode{dims_basis[1], dims_basis[2]};
    sizt_vector dims_mode_cfs{dims_weights[0]};

/*
 * Array declaration
 * --------------------------------------------
 * Name         Type                Description
 * --------------------------------------------
 * mode         Array<precision>    Individual mode, see 'lib_array.h' for datatype.
 * mode_cfs     Array<precision>    Mode-amplitudes, see 'lib_array.h' for datatype.
 */

    Array<precision> mode(dims_mode);
    Array<precision> mode_cfs(dims_mode_cfs);

/* ------------------------------------------------
 * Determine mode limit, and loop over basis modes.
 * ------------------------------------------------
 */
    sizt      mode_limit = std::min(dims_weights[0], dims_basis[0]);
    precision mode_cfs_weighted = 0.0;

    for(sizt ind = 0; ind < mode_limit; ind++){

        mode = basis.get_slice(ind, false);
    
    /* ----------------------------------------------
     * Get the mode-amplitude of the individual mode.
     * ----------------------------------------------
     */
     
        mode_cfs(ind) = (phase * mode).get_total() * (weights(ind) / norm[ind]);

    /* ------------------------------------------
     * Subtract the weighted mode from the phase.
     * ------------------------------------------
     */

        phase -= (mode * mode_cfs(ind));

    }
}

void make_psf_from_phase_screen(Array<precision>& phase, Array<precision>& psf, Array<precision>& aperture, fftw_plan& forward){

/*
 * Vector declaration
 * --------------------------------------------
 * Name             Type            Description
 * --------------------------------------------
 * dims_psf         sizt_vector     Dimensions of the PSF.
 * dims_phase       sizt_vector     Dimensions of the phase-screen.
 * dims_aperture    sizt_vector     Dimensions of the aperture.
 */

    sizt_vector dims_psf   = psf.get_dims();
    sizt_vector dims_phase = phase.get_dims();

/*
 * Variable declaration
 * --------------------------------------------
 * Name             Type            Description
 * --------------------------------------------
 * phase_center_x   sizt_vector     Center of the phase-screen.
 * phase_center_y   sizt_vector     Center of the phase-screen.
 * psf_center_x     sizt_vector     Center of the PSF.
 * psf_center_y     sizt_vector     Center of the PSF.
 */
    
    sizt phase_center_x = dims_phase[0] / 2;
    sizt phase_center_y = dims_phase[1] / 2;
    sizt psf_center_x   = dims_psf[0] / 2;
    sizt psf_center_y   = dims_psf[1] / 2;

/*
 * Variable declaration:
 * --------------------------------------------
 * Name             Type            Description
 * --------------------------------------------
 * dims_pad_begin   sizt_vector     
 */

    sizt_vector dims_pad_begin{(dims_phase[0] - dims_psf[0])/2, (dims_phase[1] - dims_psf[1])/2};

/*
 * Array declaration:
 * --------------------------------------------
 * Name             Type            Description
 * --------------------------------------------
 * pupil_function   Array<cmpx>     Array storing the pupil function.
 * complex_psf      Array<cmpx>     Array storing the complex PSF.
 *
 * --------------------
 * Additional comments:
 * --------------------
 * Pupil function has the same dimensions as the PSF.
 */

    Array<cmpx> pupil_function(dims_psf);
    Array<cmpx> complex_psf(dims_psf);

/* -------------------------
 * Construct pupil function.
 * -------------------------
 */

    for(sizt xpix = 0; xpix < dims_phase[0]; xpix++){
        for(sizt ypix = 0; ypix < dims_phase[1]; ypix++){

        /* --------------
         * Set real part.
         * --------------
         */
           
            pupil_function(xpix + phase_center_x, ypix + phase_center_y).real(static_cast<double>(aperture(xpix, ypix) * cos(phase(xpix, ypix))));
              
        /* -------------------
         * Set imaginary part.
         * -------------------
         */
           
            pupil_function(xpix + phase_center_x, ypix + phase_center_y).imag(static_cast<double>(aperture(xpix, ypix) * sin(phase(xpix, ypix)))); 
        
        }
    }

/* ------------------------------
 * Execute the fourier transform.
 * ------------------------------
 */

    fftw_execute_dft(forward, reinterpret_cast<fftw_complex*>(pupil_function[0]), reinterpret_cast<fftw_complex*>(complex_psf[0]));
    
/*
 * Array declaration.
 * --------------------------------------------
 * Name         Type                Description
 * --------------------------------------------
 * psf_copy     Array<precision>    Shifted copy of the PSF, see 'lib_array.h' for datatype.
 */

    Array<precision> psf_copy = complex_psf.get_norm().cast_to_type<precision>();

/*
 * Variable declaration.
 * --------------------------------
 * Name         Type        Description
 * --------------------------------
 * psf_total    precision   Total sum of the PSF.
 */

    precision psf_total = psf_copy.get_total();

/* ------------------------------------------------------------
 * Shift PSF central maxima to the array center, and normalize.
 * ------------------------------------------------------------
 */

    for(sizt xpix = 0; xpix < dims_psf[0]; xpix++){
        for(sizt ypix = 0; ypix < dims_psf[1]; ypix++){
            psf(xpix, ypix) = psf_copy((xpix + dims_psf[0] - psf_center_x) % dims_psf[0], (ypix + dims_psf[1] - psf_center_y) % dims_psf[1]) / psf_total;
        }
    }

}

#endif
