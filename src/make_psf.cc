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

#define _GET_ALL_POINT_SPREAD_FUNCTIONS_

/* ------------
 * Description:
 * ------------
 * This program computes the Point Spread Function(s) (PSF) of the residual phase-screens.
 * The PSF is the forward fourier transform of the pupil function, defined as:
 *
 *      Pupil_function(x, y) = aperture(x, y) * exp(i * phase(x, y))
 *
 * where 'phase' denotes an individual phase-screen.
 *
 * ------
 * Usage:
 * ------
 * mpiexec -np <cores> ./make_psf <config_file>
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
 * 2. Read phase-screen residuals from file.
 * 3. Read aperture function from file.
 * 4. Distribute residual phase-screens to MPI processes.
 * 5. Store PSF(s) returned by MPI processes.
 * 6. Repeat steps 4-5 for all residual phase-screens.
 * 7. Save PSF(s) to disk.
 *
 * -----------------------
 * Additional information:
 * -----------------------
 * Each comment block has a title that explains the succeeding code. 
 * Titles starting with !, followed by a number n indicate a block handling step n.
 */

int main(int argc, char *argv[]){

/*
 * Variable declaration:
 * --------------------------------------------
 * Name		        Type            Description
 * --------------------------------------------
 * mpi_status       MPI_status      See MPI documentation.
 * mpi_precision    MPI_Datatype    MPI_FLOAT or MPI_DOUBLE.
 */
   
    MPI_Status   mpi_status;
    MPI_Datatype mpi_precision = std::is_same<precision, float>::value == true ? MPI_FLOAT : MPI_DOUBLE;

/*
 * Variable declaration:
 * --------------------------------------------
 * Name		            Type        Description
 * --------------------------------------------
 * mpi_process_rank     int         Rank of MPI processes.
 * mpi_process_size     int         Store the total number of MPI processes
 * mpi_process_kill     int         Number of killed processes.
 * mpi_recv_count       int         Store the count of data received in MPI_Recv, see MPI documentation for explanation.
 * rd_status            int         File read status.
 * wr_status            int         File write status.
 */

    int mpi_process_rank = 0;
    int mpi_process_size = 0;
    int mpi_process_kill = 0;
    int mpi_recv_count   = 0;
    int rd_status = 0;
    int wr_status = 0;

/* --------------
 * Initialize MPI
 * --------------
 */
   
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_process_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_process_rank);

/* -----------------------------------------------------------
 * Only the master MPI process - rank zero - prints to stdout.
 * -----------------------------------------------------------
 */

    FILE *console   = mpi_process_rank == 0 ? stdout : fopen("/dev/null","wb");
    fprintf(console, "----------------------------------------------\n");
    fprintf(console, "- Point Spread Functions computation program -\n");
    fprintf(console, "----------------------------------------------\n");

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
    }
    fprintf(console, "[Done] (%s)\n", argv[1]);
    fflush (console);

/* 
 * ------------------------------
 * Workflow for root MPI process.
 * ------------------------------
 */

    if(mpi_process_rank == 0){

    /* 
     * Variable declaration.
     * ----------------------------------------
     * Name                 Type    Description
     * ----------------------------------------
     * fried_next           sizt    Index of the next fried parameter.
     * fried_done           sizt    Number of fried parameters processed.
     * percent_assigned     float   Percentage of fried assigned.
     * percent_completed    float   Percentage of fried completed.
     */

        sizt  fried_next        = 0;
        sizt  fried_done        = 0;
        float percent_assigned  = 0;
        float percent_completed = 0;

    /*
     * Array declaration.
     * --------------------------------------------
     * Name         Type                Description
     * --------------------------------------------
     * residual     Array<precision>    Phase-screen residuals, see 'lib_array.h' for datatype.
     * aperture     Array<precision>    Aperture function, see 'lib_array.h' for datatype.
     */

        Array<precision> residual;
        Array<precision> aperture;

    /* -------------------------------------------
     * !(2) Read residual phase-screens from file.
     * -------------------------------------------
     */

        fprintf(console, "(Info)\tReading file:\t\t");
        fflush (console);

        rd_status = residual.rd_fits(io_t::wr_residual_to.c_str());
        if(rd_status != EXIT_SUCCESS){
            fprintf(console, "[Failed][Err code = %d](%s)\n", rd_status, io_t::wr_residual_to.c_str());
            fflush (console);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
            
        fprintf(console, "[Done] (%s)\n", io_t::wr_residual_to.c_str());
        fflush (console);

    /* --------------------------------------
     * !(3) Read aperture function from file.
     * --------------------------------------
     */

        fprintf(console, "(Info)\tReading file:\t\t");
        fflush (console);

        rd_status = aperture.rd_fits(io_t::rd_aperture_from.c_str());
        if(rd_status != EXIT_SUCCESS){            
            fprintf(console, "[Failed][Err code = %d](%s)\n", rd_status, io_t::rd_aperture_from.c_str());
            fflush (console);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        fprintf(console, "[Done] (%s)\n", io_t::rd_aperture_from.c_str());
        fflush (console);

    /*
     * Vector declaration.
     * ------------------------------------------------
     * Name                 Type            Description
     * ------------------------------------------------
     * dims_residual        sizt_vector     Dimensions of the array storing the residual phase-screens.
     * process_fried_map    sizt_vector     Map linking process to the index of the residual phase-screen.
     */

        sizt_vector dims_residual = residual.get_dims();
        sizt_vector process_fried_map(dims_residual[0] + 1);

    /*
     * Variable declaration.
     * --------------------------------------------
     * Name             Type            Description
     * --------------------------------------------
     * dims_psf_all     sizt_vector     Dimensions of the array storing the PSFs. If the PSFs for
     *                                  all residual phase-screens is requested, then <dims_psf_all>
     *                                  represents a 4D array. If not, then it represents a 3D array.
     */

#ifdef _GET_ALL_POINT_SPREAD_FUNCTIONS_
        sizt_vector dims_psf_all{dims_residual[0], dims_residual[1], 2 * dims_residual[2] - 1, 2 * dims_residual[3] - 1};
#else
        sizt_vector dims_psf_all{dims_residual[0], 2 * dims_residual[2] - 1, 2 * dims_residual[3] - 1};
#endif

    /*
     * Array declaration.
     * --------------------------------------------
     * Name         Type                Description
     * --------------------------------------------
     * psf_all      Array<precision>   Array storing the PSFs.
     */

        Array<precision> psf_all(dims_psf_all);

    /* ----------------------------------------
     * Broadcast aperture to all MPI processes.
     * ----------------------------------------
     */
        
        MPI_Bcast(aperture[0], aperture.get_size(), mpi_precision, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

    /* -------------------------------------------------------
     * !(4) Distribute the residual phase-screens to MPI processes.
     * -------------------------------------------------------
     */

        for(int pid = 1; pid <= std::min(sizt(mpi_process_size - 1), dims_residual[0]); pid++){

        /* ---------------------------------------------------------
         * Send residual phase-screens at fried_next to MPI process.
         * ---------------------------------------------------------
         */

            MPI_Send(residual[fried_next], sizeof_vector(dims_residual, 1), mpi_precision, pid, mpi_cmds::task, MPI_COMM_WORLD);
            process_fried_map[pid] = fried_next;
            fried_next++;
        }
        percent_assigned  = (100.0 * fried_next) / dims_residual[0];
        fprintf(console, "\r(Info)\tComputing PSFs:\t\t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed); 
        fflush (console);

    /* --------------------------
     * Kill excess MPI processes.
     * --------------------------
     */        

        for(int pid = dims_residual[0] + 1; pid < mpi_process_size; pid++){   
            MPI_Send(residual[0], sizeof_vector(dims_residual, 1), mpi_precision, pid, mpi_cmds::kill, MPI_COMM_WORLD);
            mpi_process_kill++;
        }
        mpi_process_size -= mpi_process_kill;

    /* --------------------------------------------------
     * !(4) Distribute residual phase-screens to workers.
     * !(5) Store PSFs returned by workers.
     * !(6) Repeat steps 4-5 for all residual phase-screens.
     * -----------------------------------------------------
     */

        while(fried_done < dims_residual[0]){

        /* -------------------------------------------------------------------
         * Wait for an MPI process to ping root, then get process information.
         * Store the PSFs returned by the MPI process at the correct location.
         * -------------------------------------------------------------------
         */	
	
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
            MPI_Get_count(&mpi_status, mpi_precision, &mpi_recv_count);
            MPI_Recv(psf_all[process_fried_map[mpi_status.MPI_SOURCE]], sizeof_vector(dims_psf_all, 1), mpi_precision, mpi_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);

        /* --------------------------------------------------------------------------------
         * Increment <fried_done>, then display <percent_assigned> and <percent_completed>.
         * --------------------------------------------------------------------------------
         */

            fried_done++;
            percent_completed  = (100.0 * fried_done) / dims_residual[0];
            fprintf(console, "\r(Info)\tComputing PSFs:\t\t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed); 
            fflush (console);

            if(fried_next < dims_residual[0]){
        
            /* ------------------------------------------------------------------
             * If more residual phase-screens are available, send to MPI process.
             * ------------------------------------------------------------------
             */

                MPI_Send(residual[fried_next], sizeof_vector(dims_residual, 1), mpi_precision, mpi_status.MPI_SOURCE, mpi_cmds::task, MPI_COMM_WORLD);
                process_fried_map[mpi_status.MPI_SOURCE] = fried_next;

            /* --------------------------------------------------------------------------------
             * Display <percent_assigned> and <percent_completed>, then increment <fried_next>.
             * --------------------------------------------------------------------------------
             */

                percent_assigned = (100.0 * (fried_next + 1)) / dims_residual[0];
                fprintf(console, "\r(Info)\tComputing PSFs:\t\t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed); 
                fflush (console);
                fried_next++;    
            }
        }

    /* ---------------------------
     * Shutdown all MPI processes.
     * ---------------------------
     */

        for(int pid = 1; pid < mpi_process_size; pid++)
            MPI_Send(residual[0], sizeof_vector(dims_residual, 1), mpi_precision, pid, mpi_cmds::kill, MPI_COMM_WORLD);

        if(io_t::save){
    
        /* -----------------------
         * !(7) Save PSFs to disk.
         * -----------------------
         */
            
            fprintf(console, "\n(Info)\tWriting to file:\t");
            fflush (console);

            switch(format_t::wr_phase){
                case fmt_t::BIN  : wr_status = psf_all.wr_bin(io_t::wr_psf_to.c_str(), io_t::clobber);
                                   break;
                case fmt_t::FITS : wr_status = psf_all.wr_fits(io_t::wr_psf_to.c_str(), io_t::clobber);
                                   break;
                default          : wr_status = EXIT_FAILURE;
                                   break;
            }           
            if(wr_status != EXIT_SUCCESS){
                fprintf(console, "[Failed][Err code = %d](%s)\n", wr_status, io_t::wr_psf_to.c_str());
                fflush (console);
            }else{
                fprintf(console, "[Done] (%s)\n", io_t::wr_psf_to.c_str());
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
     * Vector declaration.
     * --------------------------------------------
     * Name             Type            Description
     * --------------------------------------------
     * dims_aperture    sizt_vector     Dimensions of the array storing the aperture function.
     * dims_residual    sizt_vector     Dimensions of the array storing the residual phase-screens.
     * dims_psf         sizt_vector     Dimensions of the array storing a single PSF.
     */

        const sizt_vector dims_aperture{sims_t::size_x_in_pixels, sims_t::size_y_in_pixels};
        const sizt_vector dims_residual{sims_t::realizations_per_fried, sims_t::size_x_in_pixels, sims_t::size_y_in_pixels};
        const sizt_vector dims_psf{2 * sims_t::size_x_in_pixels - 1, 2 * sims_t::size_y_in_pixels - 1}; 

    /*
     * Vector declaration:
     * --------------------------------------------
     * Name             Type            Description
     * --------------------------------------------
     * dims_psf_all     sizt_vector     Dimensions of the array storing all PSFs.
     */
       
#ifdef _GET_ALL_POINT_SPREAD_FUNCTIONS_
        const sizt_vector dims_psf_all{sims_t::realizations_per_fried, 2 * sims_t::size_x_in_pixels - 1, 2 * sims_t::size_y_in_pixels - 1};
#else
        const sizt_vector dims_psf_all(dims_psf);
#endif

    /*
     * Array declaration:
     * --------------------------------------------
     * Name         Type                Description
     * --------------------------------------------
     * aperture     Array<precision>    Aperture function.
     */

        Array<precision> aperture{dims_aperture};

    /* --------------------------------
     * Get aperture function from root.
     * --------------------------------
     */

        MPI_Bcast(aperture[0], aperture.get_size(), mpi_precision, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Array declaration.
     * --------------------------------------------------------
     * Name                     Type                Description
     * --------------------------------------------------------
     * psf                      Array<precision>    Array storing a single residual phase-screen.
     * psf_all                  Array<precision>    Array storing the residual phase-screens, per fried parameter.
     * residual                 Array<precision>    Array storing a single residual phase-screen.
     * residual_all             Array<precision>    Array storing the residual phase-screens.
     * pupil_function           Array<cmpx>         Array storing a single pupil function.
     * pupil_function_fourier   Array<cmpx>         Array storing the fourier of a single pupil function.
     */

        Array<precision> psf(dims_psf);
        Array<precision> psf_all(dims_psf_all);
        Array<precision> residual(dims_aperture);
        Array<precision> residual_all(dims_residual);
        Array<cmpx>      pupil_function(dims_psf);
        Array<cmpx>      pupil_function_fourier(dims_psf);

    /* -------------------------------
     * Import fft wisdom if available.
     * -------------------------------
     */

        fftw_import_wisdom_from_filename(io_t::rd_psf_wisdom_from.c_str());

    /*
     * Variable declaration:
     * --------------------------------
     * Name     Type        Description
     * --------------------------------
     * forward  fftw_plan   Re-usable FFTW plan for the forward transformation.
     */

        fftw_plan forward = fftw_plan_dft_2d(dims_psf[0], dims_psf[1],\
                                             reinterpret_cast<fftw_complex*>(pupil_function[0]),\
                                             reinterpret_cast<fftw_complex*>(pupil_function_fourier[0]),\
                                             FFTW_FORWARD, FFTW_ESTIMATE);

    /* ----------------------------------------------------
     * Compute PSFs of residual phase-screens until killed.
     * ----------------------------------------------------
     */

        MPI_Recv(residual_all[0], residual_all.get_size(), mpi_precision, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
            
        while(mpi_status.MPI_TAG != mpi_cmds::kill){

#ifdef _GET_ALL_POINT_SPREAD_FUNCTIONS_
            for(sizt ind = 0; ind < sims_t::realizations_per_fried; ind++){                    
                psf      = psf_all.get_slice(ind, false);
                if(!aperture_t::make_airy_disk)
                    residual = residual_all.get_slice(ind, false);
                    
                make_psf_from_phase_screen(residual, psf, aperture, forward);
            }
#else
            for(sizt ind = 0; ind < sims_t::realizations_per_fried; ind++){        
                if(!aperture_t::make_airy_disk)
                    residual = residual_all.get_slice(ind, false);
                    
                make_psf_from_phase_screen(residual, psf, aperture, forward);
                psf_all += psf;
            }
            psf_all /= sims_t::realizations_per_fried;
#endif
        /* ---------------------------------------------------------------
         * Send PSFs to MPI root, then receive new residual phase-screens.
         * ---------------------------------------------------------------
         */

            MPI_Send(psf_all[0], psf_all.get_size(), mpi_precision, 0, mpi_pmsg::ready, MPI_COMM_WORLD);
            MPI_Recv(residual_all[0], residual_all.get_size(), mpi_precision, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
        }
        
    /* -------------------------
     * Write FFT wisdom to file.
     * -------------------------
     */
   
        fftw_export_wisdom_to_filename(io_t::rd_psf_wisdom_from.c_str());
        fftw_destroy_plan(forward);
        fftw_cleanup();

    }

    MPI_Finalize();
    return(EXIT_SUCCESS);
}
