#include "mpi.h"
#include "fftw3.h"
#include "config.h"
#include "lib_mpi.h"
#include "lib_array.h"
#include "lib_phase.h"

#include <ctime>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <iostream>

/* ------------
 * Description:
 * ------------
 * This program convolves an input image with the Point Spread Functions (PSFs) corresponding to residual phase-screens.
 * The convolution is implemented as a multiplication in fourier space.
 *
 * ------
 * Usage:
 * ------
 * mpiexec -np <cores> ./make_img <config_file>
 *
 * -------
 * Inputs:
 * -------
 * See config.h for a detailed explanation of the inputs to the program.
 *
 * --------
 * Outputs:
 * --------
 * See config.h for a detailed explanation of the output of the program.
 *
 * --------------
 * Program logic:
 * -------------- 
 * 1. Read and parse the config file.
 * 2. Read the PSFs from file.
 * 3. Read the input image from file.
 * 4. Distribute PSFs to workers.
 * 5. Store the convolved images returned by workers.
 * 6. Repeat steps 4-5 for all PSFs.
 * 7. Save convolved images to disk.
 *
 * -----------------------
 * Additional information:
 * -----------------------
 * Each comment block has a title that explains the succeeding code. 
 * Titles starting with !, followed by a number n indicate a block handling step n.
 */

int main(int argc, char *argv[]){

/*
 *  Variable declaration:
 *  -------------------------------------------
 *  Name		    Type            Description
 *  -------------------------------------------
 *  mpi_status          MPI_status      See MPI documentation.
 *  mpi_precision       MPI_Datatype    MPI_FLOAT or MPI_DOUBLE.
 */
   
    MPI_Status   mpi_status;
    MPI_Datatype mpi_precision = std::is_same<precision, float>::value == true ? MPI_FLOAT : MPI_DOUBLE;

/*
 *  Variable declaration:
 *  -------------------------------------------
 *  Name		        Type        Description
 *  -------------------------------------------
 *  mpi_process_rank    int         MPI process rank.
 *  mpi_process_size    int         MPI processes, total.
 *  mpi_process_kill    int         MPI processes killed.
 *  mpi_recv_count      int         Count of data received in MPI_Recv().
 *  rd_status           int         File read status.
 *  wr_status           int         File write status.
 */

    int mpi_process_rank = 0;
    int mpi_process_size = 0;
    int mpi_process_kill = 0;
    int mpi_recv_count = 0;
    int rd_status = 0;
    int wr_status = 0;

/* ---------------
 * Initialize MPI.
 * ---------------
 */
   
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_process_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_process_rank);

/* -----------------------------------------------------------
 * Only the master MPI process - rank zero - prints to stdout.
 * -----------------------------------------------------------
 */

    FILE *console   = mpi_process_rank == 0 ? stdout : fopen("/dev/null","wb");
    fprintf(console, "-----------------------------\n");
    fprintf(console, "- Image convolution program -\n");
    fprintf(console, "-----------------------------\n");

/* ------------------------------------
 * !(1) Read and parse the config file.
 * ------------------------------------
 */

    if(argc < 2){
        fprintf(console, "(Error)\tExpected configuration file, aborting!\n");
        fflush (console);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    fprintf(console, "(Info)\tReading configuration:\t");
    fflush (console);
    
    if(config_parse(argv[1]) == EXIT_FAILURE){
        fprintf(console, "[Failed] (%s)\n", argv[1]);
        fflush (console);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }else{
        fprintf(console, "[Done] (%s)\n", argv[1]);
        fflush (console);
    }

/* 
 * ------------------
 * MPI root workflow.
 * ------------------
 */

    if(mpi_process_rank == 0){

    /* 
     * Variable declaration.
     * --------------------------------------------
     * Name                 Type    Description
     * --------------------------------------------
     * fried_next           sizt    Index of the next fried parameter.
     * fried_done           sizt    Number of fried parameters processed.
     * percent_assigned     float   Percentage of fried assigned.
     * percent_completed    float   Percentage of fried completed.
     */

        sizt  fried_next        = 0;
        sizt  fried_done        = 0;
        float percent_assigned  = 0.0;
        float percent_completed = 0.0;

    /*
     * Array declaration.
     * --------------------------------------------
     * Name     Type                Description
     * --------------------------------------------
     * img      Array<precision>    Image to be convolved, see 'lib_array.h' for datatype.
     * psfs     Array<precision>    PSFs of residual phase-screens, see 'lib_array.h' for datatype.
     */

        Array<precision> img;
        Array<precision> psfs;

    /* -----------------------------
     * !(2) Read the PSFs from file.
     * -----------------------------
     */

        fprintf(console, "(Info)\tReading file:\t\t");
        fflush (console);

        rd_status = psfs.rd_fits(io_t::wr_psf_to.c_str());
        if(rd_status != EXIT_SUCCESS){
            fprintf(console, "[Failed][Err code = %d](%s)\n", rd_status, io_t::wr_psf_to.c_str());
            fflush (console);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }else{
            fprintf(console, "[Done] (%s)\n", io_t::wr_psf_to.c_str());
            fflush (console);
        }

    /*
     * Vector declaration.
     * ----------------------------------------
     * Name         Type            Description
     * ----------------------------------------
     * dims_img     sizt_vector     Dimensions of the original image.
     */

        const sizt_vector dims_psfs = psfs.get_dims();

    /* ------------------------------
     * !(3) Read the image from file.
     * ------------------------------
     */

        fprintf(console, "(Info)\tReading file:\t\t");
        fflush (console);

        rd_status = img.rd_fits(io_t::rd_image_from.c_str());
        if(rd_status != EXIT_SUCCESS){
            fprintf(console, "[Failed][Err code = %d](%s)\n", rd_status, io_t::rd_image_from.c_str());
            fflush (console);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);        
        }else{
            fprintf(console, "[Done] (%s)\n", io_t::rd_image_from.c_str());
            fflush (console);
        }

    /*
     * Vector declaration.
     * ----------------------------------------
     * Name         Type            Description
     * ----------------------------------------
     * dims_img     sizt_vector     Dimensions of the original image.
     */

        sizt_vector dims_img = img.get_dims();

    /* --------------------------
     * Validate image dimensions.
     * --------------------------
     */

        if(dims_img.size() != 2){
            fprintf(console, "(Error)\tExpected 2D image, calling MPI_Abort()\n");
            fflush (console);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

    /*
     * Variable declaration
     * --------------------------------------------
     * Name                     Type    Description
     * --------------------------------------------
     * resampled_size_x_pixels  sizt    X-dimension of the image in pixels, after re-sampling.
     * resampled_size_y_pixels  sizt    Y-dimension of the image in pixels, after re-sampling.
     */

        sizt resampled_size_x_pixels = sizt(round(precision(dims_img[0] * image_t::original_sampling / image_t::degraded_sampling)));
        sizt resampled_size_y_pixels = sizt(round(precision(dims_img[1] * image_t::original_sampling / image_t::degraded_sampling)));

    /* -----------------------------
     * Re-sample array if requested.
     * -----------------------------
     */

        if(image_t::degraded_sampling > image_t::original_sampling){

            fprintf(console, "(Info)\tRe-sampling image:\t");
            fflush (console);

        /*
         * Vector declaration.
         * ------------------------------------------------
         * Name                     Type            Description
         * ------------------------------------------------
         * dims_shift               sizt_vector     
         * dims_img_resampled       sizt_vector
         * dims_shift_resampled     sizt_vector
         * dims_crop_coordinate     sizt_vector
         */

            const sizt_vector dims_shift{dims_img[0] / 2, dims_img[1] / 2};
            const sizt_vector dims_img_resampled{resampled_size_x_pixels, resampled_size_y_pixels};
            const sizt_vector dims_shift_resampled{dims_img_resampled[0] / 2, dims_img_resampled[1] / 2};
            const sizt_vector dims_crop_start{dims_shift[0] - dims_shift_resampled[0], dims_shift[1] - dims_shift_resampled[1]};

        /*
         * Array declaration.
         * --------------------------------------------------------
         * Name                     Type                Description
         * --------------------------------------------------------
         * img_cmpx                 Array<cmpx>
         * img_fourier              Array<cmpx>
         * img_cmpx_resampled       Array<cmpx>
         * img_fourier_resampled    Array<cmpx>
         * img_resampled            Array<precision>
         */

            Array<cmpx> img_cmpx = img.cast_to_type<cmpx>();
            Array<cmpx> img_fourier(dims_img);
            Array<cmpx> img_cmpx_resampled(dims_img_resampled);
            Array<cmpx> img_fourier_resampled(dims_img_resampled);

        /*
         * Vector declaration.
         * ------------------------------------------------
         * Name                 Type            Description
         * ------------------------------------------------
         * img_cmpx_forward     fftw_plan
         * img_cmpx_reverse     fftw_plan   
         */
                        
            fftw_plan img_cmpx_forward = fftw_plan_dft_2d(dims_img[0], dims_img[1],\
                                                          reinterpret_cast<fftw_complex*>(img_cmpx[0]),\
                                                          reinterpret_cast<fftw_complex*>(img_fourier[0]),\
                                                          FFTW_FORWARD, FFTW_ESTIMATE);

            fftw_plan img_cmpx_reverse = fftw_plan_dft_2d(dims_img_resampled[0], dims_img_resampled[1],\
                                                          reinterpret_cast<fftw_complex*>(img_fourier_resampled[0]),\
                                                          reinterpret_cast<fftw_complex*>(img_cmpx_resampled[0]),\
                                                          FFTW_BACKWARD, FFTW_ESTIMATE);
        
        /* -----------------------------------------------------------
         * Execute the forward fourier transform of the complex image.
         * Resample image in fourier space after shifting the zero frequency.
         * Execute the reverse fourier transform of the resampled complex image.
         * ---------------------------------------------------------------------
         */
        
            fftw_execute_dft(img_cmpx_forward, reinterpret_cast<fftw_complex*>(img_cmpx[0]), reinterpret_cast<fftw_complex*>(img_fourier[0]));
            img_fourier_resampled = img_fourier.get_roll(dims_shift, true).get_crop(dims_crop_start, dims_img_resampled).get_roll(dims_shift_resampled, false);
            fftw_execute_dft(img_cmpx_reverse, reinterpret_cast<fftw_complex*>(img_fourier_resampled[0]), reinterpret_cast<fftw_complex*>(img_cmpx_resampled[0]));

        /* ------------------------------------------
         * Obtain resampled image as the complex abs.
         * ------------------------------------------
         */

            img      = img_cmpx_resampled.get_abs().cast_to_type<precision>();
            dims_img = dims_img_resampled;

            fprintf(console, "[Done] (%lu %lu)\n", dims_img_resampled[0], dims_img_resampled[1]);
            fflush (console);

            fftw_destroy_plan(img_cmpx_forward);
            fftw_destroy_plan(img_cmpx_reverse);
        }

    /*
     * Vector declaration.
     * ----------------------------------------------------
     * Name                     Type            Description
     * ----------------------------------------------------
     * dims_imgs                sizt_vector     Dimensions of the array storing the convolved images.
     * dims_psfs_per_fried      sizt_vector     Dimensions of the array storing the PSFs, per fried.
     * process_fried_map        sizt_vector     Map linking an MPI process to the fried index it \\
     *                                          worked on.
     */

        sizt_vector dims_imgs(dims_psfs);
        sizt_vector dims_psfs_per_fried(dims_psfs.begin() + 1, dims_psfs.end());
        sizt_vector process_fried_map(dims_psfs[0] + 1);

        dims_imgs.rbegin()[0] = dims_img.rbegin()[0];
        dims_imgs.rbegin()[1] = dims_img.rbegin()[1];

    /*
     * Array declaration.
     * ----------------------------------------
     * Name     Type                Description
     * ----------------------------------------
     * imgs     Array<precision>    Convolved images, see 'lib_array.h' for datatype.
     */

        Array<precision> imgs(dims_imgs);

    /*
     * Variable declaration.
     * ------------------------------------------------
     * Name                         Type    Description
     * ------------------------------------------------
     * dims_psfs_per_fried_naxis     sizt    Number of dimensions of psf, per fried.
     */

        sizt dims_psfs_per_fried_naxis = dims_psfs_per_fried.size();

    /* ----------------------------------
     * Distribute image to MPI processes.
     * ----------------------------------
     */

        MPI_Bcast(&dims_psfs_per_fried_naxis,                           1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast( dims_psfs_per_fried.data(), dims_psfs_per_fried.size(), MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast( dims_img.data(),                       dims_img.size(), MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast( img[0],                                 img.get_size(), mpi_precision,     0, MPI_COMM_WORLD);

    /* ------------------------------------------
     * !(4) Distribute the PSFs to MPI processes.
     * ------------------------------------------
     */

        for(int pid = 1; sizt(pid) <= std::min(sizt(mpi_process_size) - 1, dims_psfs[0]); pid++){

        /* -------------------------------------------
         * Send the PSFs at fried_next to MPI process.
         * Record the sent index in <process_fried_map>.
         * Increment <fried_next>.
         * -----------------------
         */

            MPI_Send(psfs[fried_next], sizeof_vector(dims_psfs, 1), mpi_precision, pid, mpi_cmds::task, MPI_COMM_WORLD);
            process_fried_map[pid] = fried_next;
            fried_next++;
        }
        percent_assigned = (100.0 * (fried_next)) / dims_psfs[0];
        fprintf(console, "\r(Info)\tConvolving image:\t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed); 
        fflush (console);

    /* --------------------------
     * Kill excess MPI processes.
     * --------------------------
     */

        for(int pid = dims_psfs[0] + 1; pid < mpi_process_size; pid++){   
            MPI_Send(psfs[fried_next], sizeof_vector(dims_psfs, 1), mpi_precision, pid, mpi_cmds::kill, MPI_COMM_WORLD);
            mpi_process_kill++;
        }
        mpi_process_size -= mpi_process_kill;

    /* ------------------------------------
     * !(4) Distribute the PSFs to workers.
     * !(5) Store convolved images returned  by workers.
     * !(6) Repeat steps 4-5 for all PSFs.
     * -----------------------------------
     */

        while(fried_done < dims_psfs[0]){

        /* -------------------------------------------------------------------
         * Wait for an MPI process to ping root, then get process information.
         * Store the convolved images returned by MPI process at correct location.
         * Increment <fried_done>.
         * -----------------------
         */	
	
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
            MPI_Get_count(&mpi_status, mpi_precision, &mpi_recv_count);
            MPI_Recv(imgs[process_fried_map[mpi_status.MPI_SOURCE]], sizeof_vector(dims_imgs, 1), mpi_precision, mpi_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
            fried_done++;

        /* ---------------------------------------------------
         * Display <percent_assigned> and <percent_completed>.
         * ---------------------------------------------------
         */

            percent_completed  = (100.0 * fried_done) / dims_psfs[0];
            fprintf(console, "\r(Info)\tConvolving image:\t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed); 
            fflush (console);

            if(fried_next < dims_psfs[0]){

            /* ------------------------------------------------------------------
             * If more images need to be convolved, send new PSFs to MPI process.
             * Record sent index in <process_fried_map>.
             * -------------------------
             */

                MPI_Send(psfs[fried_next], sizeof_vector(dims_psfs_per_fried), mpi_precision, mpi_status.MPI_SOURCE, mpi_cmds::task, MPI_COMM_WORLD);
                process_fried_map[mpi_status.MPI_SOURCE] = fried_next;

            /* --------------------------------------------------------------------------------
             * Display <percent_assigned> and <percent_completed>, then increment <fried_next>.
             * --------------------------------------------------------------------------------
             */

                percent_assigned = (100.0 * (fried_next)) / dims_psfs[0];
                fprintf(console, "\r(Info)\tConvolving image:\t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed); 
                fflush (console);
                fried_next++;
            }
        }
        
        percent_assigned = (100.0 * (fried_next)) / dims_psfs[0];
        fprintf(console, "\r(Info)\tConvolving image:\t[%0.1lf %% assigned, %0.1lf %% completed]\n", percent_assigned, percent_completed); 
        fflush (console);

    /* ---------------------------
     * Shutdown all MPI processes.
     * ---------------------------
     */

        for(int pid = 1; pid < mpi_process_size; pid++)
             MPI_Send(psfs[0], sizeof_vector(dims_psfs, 1), mpi_precision, pid, mpi_cmds::kill, MPI_COMM_WORLD);

        if(io_t::save){
    
        /* -----------------------------------
         * !(7) Save convolved images to disk.
         * -----------------------------------
         */

            fprintf(console, "(Info)\tWriting to file:\t");
            fflush (console);

            switch(format_t::wr_image){
                case fmt_t::BIN  : wr_status = imgs.wr_bin(io_t::wr_image_to.c_str(), io_t::clobber);
                                   break;
                case fmt_t::FITS : wr_status = imgs.wr_fits(io_t::wr_image_to.c_str(), io_t::clobber);
                                   break;
                default          : wr_status = EXIT_FAILURE;
                                   break;
            }           
            if(wr_status != EXIT_SUCCESS){
                fprintf(console, "[Failed][Err code = %d](%s)\n", wr_status, io_t::wr_image_to.c_str());
                fflush (console);
            }else{
                fprintf(console, "[Done][%0.2lfGB](%s)\n", imgs.get_size() * sizeof(precision) / 1E9, io_t::wr_image_to.c_str());
                fflush (console);
            }
        }
    }
    
/* ---------------------
 * MPI process workflow.
 * ---------------------
 */
    
    else if(mpi_process_rank){

    /*
     * Variable declaration.
     * ------------------------------------
     * Name             Type    Description
     * ------------------------------------
     * dims_psfs_naxis  sizt    Dimensionality of the PSFs array.
     */

        sizt dims_psfs_naxis;

    /*
     * Vector declaration.
     * ----------------------------------------------------
     * Name         Type            Description
     * ----------------------------------------------------
     * dims_psf     sizt_vector     Dimensions of a single, 2D PSF. 
     * dims_img     sizt_vector     Dimensions of the 2D image to be convolved.
     * dims_psfs    sizt_vector     Dimensions of the PSFs, per fried.
     */

        sizt_vector dims_psf(2);
        sizt_vector dims_img(2);
        sizt_vector dims_psfs;

    /* ----------------------------------------------------------------------------------
     * Get dimensions of PSFs and image from root, then resize the corresponding vectors.
     * ----------------------------------------------------------------------------------
     */

        MPI_Bcast(&dims_psfs_naxis, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        dims_psfs.resize(dims_psfs_naxis);
        MPI_Bcast(dims_psfs.data(), dims_psfs.size(), MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(dims_img.data(), dims_img.size(), MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        
    /*
     * Vector declaration.
     * ----------------------------------------
     * Name         Type            Description
     * ----------------------------------------
     * dims_imgs    sizt_vector     Dimensions of the convolved images.
     * dims_shift   sizt_vector     Shift of the zero frequency in fourier array.
     */

        sizt_vector dims_imgs(dims_psfs);
        sizt_vector dims_shift{dims_img[0] / 2, dims_img[1] / 2};
 
        dims_imgs.rbegin()[0] = dims_img.rbegin()[0];
        dims_imgs.rbegin()[1] = dims_img.rbegin()[1];

        dims_psf.rbegin()[0] = dims_psfs.rbegin()[0];
        dims_psf.rbegin()[1] = dims_psfs.rbegin()[1];

    /*
     * Array declaration.
     * --------------------------------------------------------
     * Name             Type                Description
     * --------------------------------------------------------
     * img              Array<precision>    Image to be convolved.
     * imgs             Array<precision>    Image convolved for each PSF.
     * psfs             Array<precision>    Point spread functions corresponding to the residual phase-screens, per fried.
     * psf_fourier      Array<cmpx>         Fourier transform of the PSF.
     * img_fourier      Array<cmpx>         Fourier transform of the image.
     * img_fourier_c    Array<cmpx>         Fourier transform of the convolved image.
     */

        Array<precision> img (dims_img);
        Array<precision> imgs(dims_imgs);
        Array<precision> psfs(dims_psfs);
        Array<cmpx>      img_fourier(dims_img);
        Array<cmpx>      psf_fourier(dims_img);
        Array<cmpx>      img_fourier_c(dims_img);

    /* --------------------
     * Get image from root.
     * --------------------
     */

        MPI_Bcast(img[0], img.get_size(), mpi_precision, 0, MPI_COMM_WORLD);

    /*
     * Array declaration.
     * ----------------------------------------
     * Name             Type            Description
     * ----------------------------------------
     * img_double       Array<double>   Double precision image.
     * psf_double       Array<double>   Double precision PSF padded to the dimensions of the image.
     * img_double_c     Array<double>   Double precision re-sampled convolved image.
     */

        Array<double> img_double = img.cast_to_type<double>();
        Array<double> psf_double(dims_img);
        Array<double> img_double_c(dims_img);

    /*
     * Variable declaration:
     * --------------------------------
     * Name         Type        Description
     * --------------------------------
     * psf_forward  fftw_plan   Re-usable FFTW plan for the forward transformation of psf.
     * img_forward  fftw_plan   Re-usable FFTW plan for the forward transformation of img.
     * img_reverse  fftw_plan   Re-usable FFTW plan for the reverse transformation of img.
     */
    
        fftw_plan psf_forward = fftw_plan_dft_r2c_2d(dims_img[0], dims_img[1], psf_double[0], reinterpret_cast<fftw_complex*>(psf_fourier[0]), FFTW_ESTIMATE);
        fftw_plan img_forward = fftw_plan_dft_r2c_2d(dims_img[0], dims_img[1], img_double[0], reinterpret_cast<fftw_complex*>(img_fourier[0]), FFTW_ESTIMATE);
        fftw_plan img_reverse = fftw_plan_dft_c2r_2d(dims_img[0], dims_img[1], reinterpret_cast<fftw_complex*>(img_fourier_c[0]), img_double_c[0], FFTW_ESTIMATE);

    /* ----------------------------------------
     * Compute the fourier of the double image.
     * ----------------------------------------
     */

        fftw_execute(img_forward);

    /*
     * Variable declaration.
     * ----------------------------------------
     * Name             Type    Description
     * ----------------------------------------
     * img_center_x     sizt    Center of the resampled image.
     * img_center_y     sizt    Center of the resampled image.
     * psf_center_x     sizt    Center of the PSF.
     * psf_center_y     sizt    Center of the PSF.
     */

        sizt img_center_x = sizt(dims_img[0] / 2);
        sizt img_center_y = sizt(dims_img[1] / 2);
        sizt psf_center_x = sizt(dims_psf[0] / 2);
        sizt psf_center_y = sizt(dims_psf[1] / 2);

    /*
     * Variable declaration:
     * --------------------------------
     * Name         Type        Description
     * --------------------------------
     * psf_forward  fftw_plan   Re-usable FFTW plan for the forward transformation of psf.
     */

        sizt_vector dims_pad_start{(dims_img[0] - dims_psf[0]) / 2, (dims_img[1] - dims_psf[1]) / 2};

    /* -------------------------------------------
     * Convolve the image with PSFs, until killed.
     * -------------------------------------------
     */

        MPI_Recv(psfs[0], psfs.get_size(), mpi_precision, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
        while(mpi_status.MPI_TAG != mpi_cmds::kill){

        /* ----------------------------------------------------
         * Use time to seed the poisson random number generator
         * ----------------------------------------------------
         */

            Time::time_point time_start = Time::now();
            if(dims_psfs.size() == 2){
                    
            /* ------------------------------------
             * Pad the PSFs to the size of the img.
             * ------------------------------------
             */

                psf_double = psfs.cast_to_type<double>().get_pad(dims_pad_start, dims_img);
                  
            /* -----------------------------------------
             * Execute the forward transform of the PSF.
             * Convolve the image in fourier space.
             * Execute the reverse transform of the convolved image.
             * -----------------------------------------------------
             */

                fftw_execute_dft_r2c(psf_forward, psf_double[0], reinterpret_cast<fftw_complex*>(psf_fourier[0]));
                img_fourier_c = img_fourier * psf_fourier;
                fftw_execute_dft_c2r(img_reverse, reinterpret_cast<fftw_complex*>(img_fourier_c[0]), img_double_c[0]);

            /* ------------------------------------------------------------------------
             * Shift the convolved image and normalize by the total number of elements.
             * The cast the convolved image to type <precision>.
             * -------------------------------------------------
             */

                img_double_c = img_double_c.get_roll(dims_shift) / img_double_c.get_size();
                imgs = img_double_c.cast_to_type<precision>();   

            }else if(dims_psfs.size() == 3){

                for(sizt ind = 0; ind < dims_psfs[0]; ind++){

                /* ------------------------------------
                 * Pad the PSFs to the size of the img.
                 * ------------------------------------
                 */

                    psf_double = psfs.get_slice(ind).cast_to_type<double>().get_pad(dims_pad_start, dims_img);

                /* -----------------------------------------
                 * Execute the forward transform of the PSF.
                 * Convolve the image in fourier space.
                 * Execute the reverse transform of the convolved image.
                 * -----------------------------------------------------
                 */

                    fftw_execute_dft_r2c(psf_forward, psf_double[0], reinterpret_cast<fftw_complex*>(psf_fourier[0]));
                    img_fourier_c = img_fourier * psf_fourier;
                    fftw_execute_dft_c2r(img_reverse, reinterpret_cast<fftw_complex*>(img_fourier_c[0]), img_double_c[0]);

                /* -------------------------------------------------------------------------
                 * Shift the convolved image, and normalize by the total number of elements.
                 * -------------------------------------------------------------------------
                 */
                    
                    img_double_c = img_double_c.get_roll(dims_shift) / img_double.get_size();
                    for(sizt xpix = 0; xpix < dims_img[0]; xpix++)
                        for(sizt ypix = 0; ypix < dims_img[1]; ypix++)
                            imgs(ind, xpix, ypix) = static_cast<precision>(img_double_c(xpix, ypix));

                /*
                 * Variable declaration
                 * --------------------------------------------
                 * Name         Type                Description
                 * --------------------------------------------
                 * duration     Time::duration      Time duration until execution of this line. 
                 * mersenne     sizt                Seed for mersenne twister random number generator.
                 * gemerator    std::mt19937        Mersenne Twister Random number generator. 
                 */
    
                    Time::duration duration = Time::now() - time_start;
                    sizt           mersenne = duration.count();
                    std::mt19937   generator(mersenne);

                /* --------------------------------------
                 * Add the noise to the convolved images.
                 * --------------------------------------
                 */

                    for(sizt xpix = 0; xpix < dims_img[0]; xpix++){
                        for(sizt ypix = 0; ypix < dims_img[1]; ypix++){
                            std::poisson_distribution<int> distribution(image_t::normalization * imgs(ind, xpix, ypix));
                            imgs(ind, xpix, ypix) = distribution(generator);
                        }
                    }
                }
            }
            
            MPI_Send(imgs[0], imgs.get_size(), mpi_precision, 0, mpi_pmsg::ready, MPI_COMM_WORLD);
            MPI_Recv(psfs[0], psfs.get_size(), mpi_precision, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
        }

        fftw_destroy_plan(psf_forward);
        fftw_destroy_plan(img_forward);
        fftw_destroy_plan(img_reverse);
        fftw_cleanup();
    }

    MPI_Finalize();
    return(EXIT_SUCCESS);
}
