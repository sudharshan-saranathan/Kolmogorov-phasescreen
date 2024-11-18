#include "mpi.h"
#include "fftw3.h"
#include "fitsio.h"
#include "config.h"
#include "lib_mpi.h"
#include "lib_array.h"
#include "lib_phase.h"

#include <ctime>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <iostream>

#define _APERTURE_

/* ------------
 * Description:
 * ------------
 * This program subtracts basis functions defined on the aperture from phase_all-screen simulations.
 * The subtraction is weighted, if coeff are provided. The residual phase-screens represent
 * corrections to the phase-screens either by Adaptive Optics (AO) or by post-facto image
 * processing. 
 *
 * ------
 * Usage:
 * ------
 * mpiexec -np <cores> ./make_residual <config_file>
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
 * 2. Read phase_all-screen simulations from file.
 * 3. Read basis functions from file.
 * 4. Read basis coeff from file.
 * 5. Distribute phase_all-screen simulations to workers.
 * 6. Store residual phase-screens returned by workers.
 * 7. Repeat steps 5-6 for all phase_all-screen simulations.
 * 8. Save residual phase-screens to disk.
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
 *  Name            Type            Description
 *  -------------------------------------------
 *  mpi_status      MPI_status      See MPI documentation.
 *  mpi_precision   MPI_Datatype    MPI_FLOAT or MPI_DOUBLE.
 */
   
    MPI_Status   mpi_status;
    MPI_Datatype mpi_precision = std::is_same<precision, float>::value == true ? MPI_FLOAT : MPI_DOUBLE;

/*
 * Variable declaration:
 * ----------------------------------------
 * Name                 Type    Description
 * ----------------------------------------
 * mpi_process_rank     int     Rank of MPI process, see MPI documentation.
 * mpi_process_size     int     Number of MPI processes, total.
 * mpi_process_kill     int     Number of MPI processes, killed.
 * mpi_recv_count       int     Count of data received in MPI_Recv(), see MPI documentation.
 * rd_status            int     Read  status of file.
 * wr_status            int     Write status of file.
 */

    int mpi_process_rank = 0;
    int mpi_process_size = 0;
    int mpi_recv_count   = 0;
    int rd_status        = 0;
    int wr_status        = 0;

/* -------------------------
 * Initialize MPI framework.
 * -------------------------
 */
   
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_process_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_process_rank);

/* ----------------------------------------------------
 * Only the root MPI process (rank 0) prints to stdout.
 * ----------------------------------------------------
 */

    FILE *console   = mpi_process_rank == 0 ? stdout : fopen("/dev/null","wb");
    fprintf(console, "------------------------------------------------------\n");
    fprintf(console, "- phase_all-screen residuals computation program -\n");
    fprintf(console, "------------------------------------------------------\n");

/* ------------------------------------
 * !(1) Read and parse the config file.
 * ------------------------------------
 */

    if(argc < 2){
        fprintf(console, "(Error)\tExpected configuration file, calling MPI_Abort()\n");
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
 * ------------------
 * MPI root workflow.
 * ------------------
 */

    if(mpi_process_rank == 0){

    /*
     * Variable declaration:
     * ----------------------------------------
     * Name                 Type    Description
     * ----------------------------------------
     * fried_next           sizt    Index of the next fried parameter.
     * fried_done           sizt    Number of fried parameters processed.
     * percent_assigned     float   Percentage of fried parameters assigned to MPI processs.
     * percent_completed    float   Percentage of fried parameters completed by MPI processs.
     */

        sizt  fried_next        = 0;
        sizt  fried_done        = 0;
        float percent_assigned  = 0.0;
        float percent_completed = 0.0;

    /*
     * Array declaration.
     * --------------------------------------------
     * Name         Type                Description
     * ---------------------------------------------
     * basis        Array<precision>    Array storing the basis functions defined on the aperture.
     * coeff      Array<precision>    Array storing the coeff to subtract the basis functions with.
     * phase_all    Array<precision>    Array storing all phase_all-screen simulations.
     */

        Array<precision> basis;
        Array<precision> coeff;
        Array<precision> phase_all;

    /* ---------------------------------------------
     * !(2) Read phase-screen simulations from file.
     * ---------------------------------------------
     */

        fprintf(console, "(Info)\tReading file:\t\t");
        fflush (console);

        rd_status = phase_all.rd_fits(io_t::wr_phase_to.c_str());
        if(rd_status != EXIT_SUCCESS){
            fprintf(console, "[Failed][Err code = %d](%s)\n", rd_status, io_t::wr_phase_to.c_str());
            fflush (console);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        fprintf(console, "[Done] (%s)\n", io_t::wr_phase_to.c_str());
        fflush (console);

    /* ------------------------------------
     * !(3) Read basis functions from file.
     * ------------------------------------
     */

        fprintf(console, "(Info)\tReading file:\t\t");
        fflush (console);

        rd_status = basis.rd_fits(io_t::rd_basis_from.c_str());
        if(rd_status != EXIT_SUCCESS){        
            fprintf(console, "[Failed][Err code = %d](%s)\n", rd_status, io_t::rd_basis_from.c_str());
            fflush (console);        
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        fprintf(console, "[Done] (%s)\n", io_t::rd_basis_from.c_str());
        fflush (console);

    /*
     * Vector declaration.
     * ------------------------------------------------
     * Name                 Type            Description
     * ------------------------------------------------
     * dims_basis           sizt_vector     Dimensions of the array storing the basis functions.
     * dims_coeff         sizt_vector     Dimensions of the array storing the basis coeff.
     * dims_phase_all       sizt_vector     Dimensions of the array storing the phase-screens.
     * process_fried_map    sizt_vector     Map linking an MPI process to the index of fried parameter.
     */

        sizt_vector dims_basis       = basis.get_dims();
        sizt_vector dims_phase_all   = phase_all.get_dims();
        sizt_vector dims_coeff{dims_phase_all[0], dims_basis[0]};
        sizt_vector process_fried_map(dims_phase_all[0] + 1);

    /* -----------------------------------------------
     * !(4) Read basis coeff from file, if provided.
     * -----------------------------------------------
     */
        if(io_t::rd_coeff_from == "_ONES_"){

            fprintf(console, "(Info)\tUsing coeff = 1\t");
            fflush (console);

        /* -----------------------
         * Set all coeff to 1.0.
         * -----------------------
         */
            Array<precision> AO_coeff(dims_coeff);
            for(sizt xpix = 0; xpix < dims_coeff[0]; xpix++)
                for(sizt ypix = 0; ypix < dims_coeff[1]; ypix++)
                    AO_coeff(xpix, ypix) = 1.0;

            coeff = AO_coeff;

        }else if(io_t::rd_coeff_from == "_ZEROS_"){

            fprintf(console, "(Info)\tUsing coeff = 0\t");
            fflush (console);

        /* -------------------------------------------------------------
         * Set all coeff to 0.0 (automatic when creating a new array).
         * -------------------------------------------------------------
         */
            
            Array<precision> AO_coeff(dims_coeff);
            coeff = AO_coeff;

        }else{

            fprintf(console, "(Info)\tReading file:\t\t");
            fflush (console);
         
            rd_status = coeff.rd_fits(io_t::rd_coeff_from.c_str());
            if(rd_status != EXIT_SUCCESS){        
                fprintf(console, "[Failed][Err code = %d](%s)\n", rd_status, io_t::rd_coeff_from.c_str());
                fflush (console);                
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);        
            }else if(dims_coeff[0] != coeff.get_dims(0)){
                fprintf(console, "(Error)\tExpected coeff with dimensions [%lu, %lu]\n", dims_coeff[0], dims_coeff[1]);
                fflush (console);
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }            
            dims_coeff = coeff.get_dims();
            fprintf(console, "[Done] (%s)\n", io_t::rd_coeff_from.c_str());
            fflush (console);
        }

    /* --------------------------------------------
     * Validate the dimensions of the input arrays.
     * --------------------------------------------
     */

        if(dims_phase_all[2] != sims_t::size_x_in_pixels  || dims_phase_all[3] != sims_t::size_y_in_pixels){     
            fprintf(console, "(Error)\tExpected phase-screens with size [%lu %lu], calling MPI_Abort()\n", sims_t::size_x_in_pixels, sims_t::size_y_in_pixels);
            fflush (console);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        if(dims_basis[1] != sims_t::size_x_in_pixels  || dims_basis[2] != sims_t::size_y_in_pixels){     
            fprintf(console, "(Error)\tExpected basis functions with size [%lu %lu], calling MPI_Abort()\n", sims_t::size_x_in_pixels, sims_t::size_y_in_pixels);
            fflush (console);            
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

     /*
     * Array declaration:
     * --------------------------------------------
     * Name         Type                Description
     * --------------------------------------------
     * residual     Array<precision>    Residual phase-screens, see 'lib_array.h' for datatype.
     */

        Array<precision> residual(phase_all);

    /* -----------------------------------------------
     * Broadcast basis functions to all MPI processes.
     * -----------------------------------------------
     */

        MPI_Bcast(dims_basis.data(), dims_basis.size(), MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(basis[0], basis.get_size(), mpi_precision, 0, MPI_COMM_WORLD);

    /* --------------------------------------------------------
     * !(5) Distribute phase_all-screen simulations to workers.
     * !(6) Store residual phase-screens returned by workers.
     * !(7) Repeat steps 5-6 for all phase_all-screen simulations.
     * -----------------------------------------------------------
     */

        for(int pid = 1; pid <= std::min(sizt(mpi_process_size - 1), dims_phase_all[0]); pid++){

        /* -------------------------------------------------------------------
         * Send phase-screens simulations, and basis coeff to MPI processes.
         * Record the index sent in <process_fried_map>.
         * ---------------------------------------------
         */

            MPI_Send(phase_all[fried_next], sizeof_vector(dims_phase_all, 1), mpi_precision, pid, mpi_cmds::task, MPI_COMM_WORLD);
            MPI_Send(coeff[fried_next], dims_coeff[1], mpi_precision, pid, mpi_cmds::task, MPI_COMM_WORLD);
            process_fried_map[pid] = fried_next;
            fried_next++;
        }
        percent_assigned  = (100.0 * fried_next) / dims_phase_all[0];
        fprintf(console, "\r(Info)\tComputing residuals: \t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed); 
        fflush (console);

    /* --------------------------
     * Kill excess MPI processes.
     * --------------------------
     */

        for(int pid = dims_phase_all[0] + 1; pid < mpi_process_size; pid++){
            MPI_Send(phase_all[0], sizeof_vector(dims_phase_all, 1), mpi_precision, pid, mpi_cmds::kill, MPI_COMM_WORLD);
            MPI_Send(coeff[0], dims_coeff[1], mpi_precision, pid, mpi_cmds::kill, MPI_COMM_WORLD);
        }

        while(fried_done < dims_phase_all[0]){

        /* ---------------------------------------------------------------------
         * Wait for an MPI process to ping master, then get process information.
         * Store the residual phase-screens returned by MPI process at the correct location.
         * Increment <fried_done>.
         * -----------------------
         */	
	
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
            MPI_Get_count(&mpi_status, mpi_precision, &mpi_recv_count);
            MPI_Recv(residual[process_fried_map[mpi_status.MPI_SOURCE]], sizeof_vector(dims_phase_all, 1), mpi_precision, mpi_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
            fried_done++;

        /* ---------------------------------------------------
         * Display <percent_assigned> and <percent_completed>.
         * ---------------------------------------------------
         */

            percent_completed  = (100.0 * fried_done) / dims_phase_all[0];
            fprintf(console, "\r(Info)\tComputing residuals:\t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed); 
            fflush (console);


            if(fried_next < dims_phase_all[0]){

           /* ---------------------------------------------------------------------------------
            * If more residuals need to be computed, send next set of phase-screen simulations.
            * Send the corresponding coeff, and record the index sent in <process_fried_map>.
            * ---------------------------------------------------------------------------------
            */

                MPI_Send(phase_all[fried_next], sizeof_vector(dims_phase_all, 1), mpi_precision, mpi_status.MPI_SOURCE, mpi_cmds::task, MPI_COMM_WORLD);
                MPI_Send(coeff[fried_next], dims_coeff[1], mpi_precision, mpi_status.MPI_SOURCE, mpi_cmds::task, MPI_COMM_WORLD);
                process_fried_map[mpi_status.MPI_SOURCE] = fried_next;

            /* --------------------------------------------------------------------------------
             * Display <percent_assigned> and <percent_completed>, then increment <fried_next>.
             * --------------------------------------------------------------------------------
             */

                percent_assigned = (100.0 * (fried_next + 1)) / dims_phase_all[0];
                fprintf(console, "\r(Info)\tComputing residuals:\t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed); 
                fflush (console);
                fried_next++;
            }    
        }

    /* ---------------------------
     * Shutdown all MPI processes.
     * ---------------------------
     */

        for(int pid = 1; pid < mpi_process_size; pid++){
            MPI_Send(phase_all[0], sizeof_vector(dims_phase_all, 1), mpi_precision, pid, mpi_cmds::kill, MPI_COMM_WORLD);
            MPI_Send(coeff[0], dims_coeff[1], mpi_precision, pid, mpi_cmds::kill, MPI_COMM_WORLD);
        }

        if(io_t::save){
    
        /* -----------------------------------------
         * !(8) Save residual phase-screens to disk.
         * -----------------------------------------
         */

            fprintf(console, "\n(Info)\tWriting to file:\t");
            fflush (console);

            switch(format_t::wr_residual){
                case fmt_t::BIN  : wr_status = residual.wr_bin(io_t::wr_residual_to.c_str(), io_t::clobber);
                                   break;
                case fmt_t::FITS : wr_status = residual.wr_fits(io_t::wr_residual_to.c_str(), io_t::clobber);
                                   break;
                default          : wr_status = EXIT_FAILURE;
                                   break;
            }           
            if(wr_status != EXIT_SUCCESS){
                fprintf(console, "[Failed][Err code = %d](%s)\n", wr_status, io_t::wr_residual_to.c_str());
                fflush (console);
            }else{
                fprintf(console, "[Done] (%s)\n", io_t::wr_residual_to.c_str());
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
     * dims_basis       sizt_vector     Dimensions of the array storing the basis functions.
     */
        
        sizt_vector dims_basis(3);

    /* ------------------------------------------------------
     * Get dimensions of the basis functions array from root.
     * ------------------------------------------------------
     */

        MPI_Bcast(dims_basis.data(), dims_basis.size(), MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    /* ------------------------------------
     * Get basis functions array from root.
     * ------------------------------------
     */
     
        Array<precision> basis(dims_basis);
        MPI_Bcast(basis[0], basis.get_size(), mpi_precision, 0, MPI_COMM_WORLD);
    
    /*
     * Vector declaration.
     * ----------------------------------------------------
     * Name             Type            Description
     * ----------------------------------------------------
     * dims_modes       sizt_vector     Number of basisfunctions.
     * dims_phase       sizt_vector     Dimensions of the array storing a single phase-screen.
     * dims_phase_all   sizt_vector     Dimensions of the array storing all phase-screens, per fried.
     */

        sizt_vector dims_modes{dims_basis[0]};
        sizt_vector dims_phase{sims_t::size_x_in_pixels, sims_t::size_y_in_pixels};
        sizt_vector dims_phase_all{sims_t::realizations_per_fried, dims_phase[0], dims_phase[1]};

    /*
     * Array declaration.
     * --------------------------------------------
     * Name         Type                Description
     * --------------------------------------------
     * coeff      Array<precision>    Array storing the coeff of the basis functions.
     * phase        Array<precision>    Array storing a single phase-screen.
     * phase_all    Array<precision>    Array storing the phase-screen simulations.
     */

        Array<precision> coeff(dims_modes);
        Array<precision> phase(dims_phase);
        Array<precision> phase_all(dims_phase_all);

    /*
     * Vector declaration.
     * ----------------------------------------
     * Name         Type            Description
     * ----------------------------------------
     * basis_norm   sizt_vector     L2 normalization of the basis functions.
     */

        sizt_vector basis_norm(dims_basis[0]);

    /* ------------------------------------------------
     * Compute L2 normalization of the basis functions.
     * ------------------------------------------------
     */

        for(sizt ind = 0; ind < dims_basis[0]; ind++)
            for(sizt xs = 0; xs < dims_basis[1]; xs++)
                for(sizt ys = 0; ys < dims_basis[2]; ys++)
                    basis_norm[ind] += basis(ind, xs, ys) * basis(ind, xs, ys);

    /* -----------------------------------------------------------
     * Enter loop to compute residual phase-screens, until killed.
     * -----------------------------------------------------------
     */
        
        MPI_Recv(phase_all[0], phase_all.get_size(), mpi_precision, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
        MPI_Recv(coeff[0], dims_modes[0], mpi_precision, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
        while(mpi_status.MPI_TAG != mpi_cmds::kill){
            for(sizt ind = 0; ind < dims_phase_all[0]; ind++){
                phase = phase_all.get_slice(ind, false);
                make_residual_phase_screen(phase, basis, coeff, basis_norm);
            }

            /* --------------------------------------------------------------------------------
             * Send residual phase-screens to root, then receive new phase-screens and coeff.
             * --------------------------------------------------------------------------------
             */

            MPI_Send(phase_all[0], phase_all.get_size(), mpi_precision, 0, mpi_pmsg::ready, MPI_COMM_WORLD);
            MPI_Recv(phase_all[0], phase_all.get_size(), mpi_precision, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
            MPI_Recv(coeff[0], dims_modes[0], mpi_precision, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
        }
    }

    MPI_Finalize();
    return(EXIT_SUCCESS);
}
