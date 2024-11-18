#include "mpi.h"
#include "fftw3.h"
#include "config.h"
#include "fitsio.h"
#include "lib_mpi.h"
#include "lib_array.h"
#include "lib_phase.h"

#include <ctime>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <unistd.h>
#include <iostream>

/* ------------
 * Description:
 * ------------
 *
 * This program simulates the degradation of phase-screens by atmospheric turbulence. The   \\
 * simulated phase-screens, therefore, statistically represent the spatial distribution of  \\
 * wave-front errors at the exit pupil of a telescope. The phase-screens are simulated with \\
 * the property that their power spectrum follows the Kolmogorov-Obukhov power law. Each    \\
 * phase-screen realization is computed as the fourier transform of a fourier array.
 *
 * ------
 * Usage:
 * ------
 * mpiexec -np <cores> ./make_phase <config_file>
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
 * 2. Read Fried parameters from file.
 * 3. Distribute the Fried parameters to MPI processs.
 * 4. Store the simulated phase-screens returned by MPI processs.
 * 5. Repeat steps 3-4 for all fried parameters.
 * 6. Save simulations to disk.
 *
 * -----------------------
 * Additional information:
 * -----------------------
 * Each comment block has a title that explains the succeeding code. 
 * Titles starting with !, followed by a number n indicate a block handling step n.
 */

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::duration<float>       Period;

int main(int argc, char *argv[]){

/*
 * Variable declaration:
 * --------------------------------------------
 * Name             Type            Description
 * --------------------------------------------
 * mpi_status       MPI_status      MPI status, see MPI documentation.
 * mpi_precision    MPI_Datatype    MPI Datatype, see MPI documentation.
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
    int mpi_process_kill = 0;
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

/* ------------------------------------------------------
 * Only the master MPI process (rank 0) prints to stdout.
 * ------------------------------------------------------
 */

    FILE *console   = mpi_process_rank == 0 ? stdout : fopen("/dev/null","wb");
    fprintf(console, "------------------------------------------------------\n");
    fprintf(console, "- Turbulence-degraded phasescreen simulation program -\n");
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

/* ------------------
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
        float percent_assigned  = 0;
        float percent_completed = 0;


    /* -------------------------------------
     * !(2) Read fried parameters from file.
     * -------------------------------------
     */

        fprintf(console, "(Info)\tReading file:\t\t");
        fflush (console);

    /*
     * Array declaration:
     * ----------------------------------------
     * Name		Type		        Description
     * ----------------------------------------
     * fried	Array<precision>	Array storing the fried parameters.
     */
 
        Array<precision> fried;
        rd_status = fried.rd_fits(io_t::rd_fried_from.c_str());
        if(rd_status != EXIT_SUCCESS){
            fprintf(console, "[Failed][Err code = %d](%s)\n", rd_status, io_t::rd_fried_from.c_str());
            fflush (console);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        fprintf(console, "[Done] (%s)\n", io_t::rd_fried_from.c_str());
        fflush (console);

    /*
    * Vector declaration:
    * -------------------------------------------------
    * Name			        Type			Description
    * -------------------------------------------------
    * dims_phase_all        sizt_vector     Dimensions of the array storing the phase-screens.
    * process_fried_map     sizt_vector     Map linking an MPI process to the index of fried parameter.
    */
    
        sizt_vector dims_phase_all{fried.get_size(), sims_t::realizations_per_fried, sims_t::size_x_in_pixels, sims_t::size_y_in_pixels};
        sizt_vector process_fried_map(fried.get_size() + 1);
    
#ifdef _USE_APERTURE_
    
    /* -----------------------------------------------
     * If aperture function available, read from file.
     * -----------------------------------------------
     */
	
        fprintf(console, "(Info)\tReading file:\t\t");
        fflush (console);

        Array<precision> aperture;
        rd_status = aperture.rd_fits(io_t::rd_aperture_from.c_str());
        if(rd_status != EXIT_SUCCESS){
            fprintf(console, "[Failed][Err code = %d](%s)\n", rd_status, io_t::rd_aperture_from.c_str());
            fflush (console);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	    }

    /* --------------------------------------------------------------------------
     * Check if dimensions of aperture match the values specified in config file.
     * --------------------------------------------------------------------------
     */

        if(aperture.get_dims(0) != sims_t::size_x_in_pixels && aperture.get_dims(1) != sims_t::size_y_in_pixels){
            fprintf(console, "[Failed][Expected aperture size = [%lu, %lu]]\n", sims_t::size_x_in_pixels, sims_t::size_y_in_pixels);
            fflush (console);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	    }
        fprintf(console, "[Done] (%s)\n", io_t::rd_aperture_from.c_str());
        fflush (console);
   
    /* ----------------------------------------
     * Broadcast aperture to all MPI processes.
     * ----------------------------------------
     */

        MPI_Bcast(aperture[0], aperture.get_size(), mpi_precision, 0, MPI_COMM_WORLD);
#endif 

    /*
     * Array declaration:
     * --------------------------------------------
     * Name         Type                Description
     * --------------------------------------------
     * phase_all	Array<precision>	Array storing all phase-screen simulations.
     */

        Array<precision> phase_all(dims_phase_all);

    /* ------------------------------------------------------
     * !(3) Distribute the fried parameters to MPI processes.
     * !(4) Store the simulated phase-screens returned by MPI processes.
     * !(5) Repeat steps 3-4 for all fried parameters.
     * -----------------------------------------------
     */

        percent_assigned  = (100.0 * fried_next) / fried.get_size();
        percent_completed = (100.0 * fried_done) / fried.get_size();
        
        fprintf(console, "\r(Info)\tSimulating phases:\t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed);
        fflush (console);

     /* --------------------------------------------------------
      * For pid <= fried parameters, distribute fried parameter.
      * --------------------------------------------------------
      */

        Time::time_point time_start = Time::now();
        for(int pid = 1; pid <= std::min(sizt(mpi_process_size) - 1, fried.get_size()); pid++){

        /* ------------------------------------
         * Send fried parameter to MPI process.
         * Record index sent in <process_fried_map>.
         * Increment <fried_next>.
         * -----------------------
         */

            MPI_Send(fried[fried_next], 1, mpi_precision, pid, mpi_cmds::task, MPI_COMM_WORLD);
            process_fried_map[pid] = fried_next;
            fried_next++;
        }
        percent_assigned  = (100.0 * fried_next) / fried.get_size();
        fprintf(console, "\r(Info)\tSimulating phases:\t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed); 
        fflush (console);

     /* ---------------------------------------------
      * For pid > fried parameters, kill MPI process.
      * ---------------------------------------------
      */

        for(int pid = fried.get_size() + 1; pid < mpi_process_size; pid++){   
            MPI_Send(fried[0], 1, mpi_precision, pid, mpi_cmds::kill, MPI_COMM_WORLD);
            mpi_process_kill++;
        }
        mpi_process_size -= mpi_process_kill;

    /* --------------------------------------------------------
     * Loop to simulate phase-screens for all fried parameters.
     * --------------------------------------------------------
     */

        while(fried_done < fried.get_size()){
        	  
        /* ------------------------------------------------------------------
         * Wait for a MPI process to ping root, then get process information.
         * Store phase-screen simulations returned by MPI process at the correct location.
         * Increment <fried_done>.
         * -------------------------------------------------------
         */	
	
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
            MPI_Get_count(&mpi_status, mpi_precision, &mpi_recv_count);
            MPI_Recv(phase_all[process_fried_map[mpi_status.MPI_SOURCE]], sizeof_vector(dims_phase_all, 1),\
                     mpi_precision, mpi_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
            fried_done++;
      
        /* ---------------------------------------------------
         * Display <percent_assigned> and <percent_completed>.
         * ---------------------------------------------------
         */

            percent_completed = (100.0 * fried_done) / fried.get_size();	        
            fprintf(console, "\r(Info)\tSimulating phases:\t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed);
            fflush (console);

            if(fried_next < fried.get_size()){
        
            /* ------------------------------------------------------------
             * If more fried parameters are available, send to MPI process.
             * ------------------------------------------------------------
             */
	    
                MPI_Send(fried[fried_next], 1, mpi_precision, mpi_status.MPI_SOURCE, mpi_cmds::task, MPI_COMM_WORLD);
                process_fried_map[mpi_status.MPI_SOURCE] = fried_next;
                fried_next++;

            /* ---------------------------------------------------
             * Display <percent_assigned> and <percent_completed>.
             * ---------------------------------------------------
             */

                percent_assigned  = (100.0 * fried_next) / fried.get_size();
                fprintf(console, "\r(Info)\tSimulating phases:\t[%0.1lf %% assigned, %0.1lf %% completed]", percent_assigned, percent_completed); 
                fflush (console); 
            }
        }

        Period duration = Time::now() - time_start;
        fprintf(console, "\r(Info)\tSimulating phases:\t[%0.1lf %% assigned, %0.1lf %% completed] (%0.2lf)", percent_assigned, percent_completed, duration.count()); 
    
    /* ---------------------------
     * Shutdown all MPI processes.
     * ---------------------------
     */

        for(int pid = 1; pid < mpi_process_size; pid++)
            MPI_Send(fried[0], 1, mpi_precision, pid, mpi_cmds::kill, MPI_COMM_WORLD);


        if(io_t::save){
    
        /* ------------------------------ 
         * !(6) Save simulations to disk.
         * ------------------------------
         */

            fprintf(console, "\n(Info)\tWriting to file:\t");
            fflush (console);

            switch(format_t::wr_phase){
                case fmt_t::BIN  : wr_status = phase_all.wr_bin(io_t::wr_phase_to.c_str(), io_t::clobber);
                                   break;
                case fmt_t::FITS : wr_status = phase_all.wr_fits(io_t::wr_phase_to.c_str(), io_t::clobber);
                                   break;
                default          : wr_status = EXIT_FAILURE;
                                   break;
            }           
            if(wr_status != EXIT_SUCCESS){
                fprintf(console, "[Failed][Err code = %d](%s)\n", wr_status, io_t::wr_phase_to.c_str());
                fflush (console);
            }else{
                fprintf(console, "[Done] (%s)\n", io_t::wr_phase_to.c_str());
                fflush (console);
            }
        }
    /*
     * -----------------------------
     * End of workflow for MPI root.
     * -----------------------------
     */
    }
    
/*
 * ---------------------
 * MPI process workflow.
 * ---------------------
 */    
    
    else if(mpi_process_rank){
    
    /*
     * Vector declaration:
     * ------------------------------------------------
     * Name                 Type            Description
     * ------------------------------------------------
     * dims_phase_all       sizt_vector     Dimensions of the array of all phase-screens, per fried.
     * dims_phase_cropped   sizt_vector     Dimensions of the array of a single cropped phase-screen.
     * dims_phase_complex   sizt_vector     Dimensions of the array of a single simulated phase-screen.
     *
     * --------------------
     * Additional comments:
     * --------------------
     * The physical size of the simulation should be much larger than the physical size of the \\
     * cropped phase-screen, in order to adequately reproduce the power in lower spatial orders.
     */

        const sizt_vector dims_phase_all{sims_t::realizations_per_fried, sims_t::size_x_in_pixels, sims_t::size_y_in_pixels};
        const sizt_vector dims_phase_cropped{sims_t::size_x_in_pixels, sims_t::size_y_in_pixels};
        const sizt_vector dims_phase_complex{sizt(sims_t::size_in_meters * sims_t::size_x_in_pixels * aperture_t::sampling_factor / aperture_t::size),\
                                             sizt(sims_t::size_in_meters * sims_t::size_y_in_pixels * aperture_t::sampling_factor / aperture_t::size)};

    /*
     * Vector declaration
     * --------------------------------------------
     * Name             Type            Description
     * --------------------------------------------
     * dims_crop_begin  sizt_vector     The starting coordinate for cropping the simulations.
     */

        const sizt_vector dims_crop_begin{(dims_phase_complex[0] - dims_phase_cropped[0]) / 2, (dims_phase_complex[1] - dims_phase_cropped[1]) / 2};

    /*
     * Array declaration:
     * ------------------------------------------------
     * Name             Type                Description
     * ------------------------------------------------
     * phase_all        Array<precision>    Array storing all phase-screens, per fried.
     * phase_cropped    Array<precision>    Array storing a single cropped phase-screen.
     * phase_complex    Array<cmpx>         Array storing a single simulated phase-screen.
     * phase_fourier	Array<cmpx>         Array storing the fourier of a single simulated phase-screen.
     *
     * --------------------
     * Additional comments:
     * --------------------
     * phase_complex and phase_fourier are re-used over the requested number of realizations, 
     * cropped to the dimensions of the aperture and stored in phase_per_fried. 
     * phase_all is then sent to MPI rank = 0.
     */

        Array<precision> phase_all(dims_phase_all);
        Array<precision> phase_cropped(dims_phase_cropped);
        Array<cmpx>      phase_complex(dims_phase_complex);
        Array<cmpx>      phase_fourier(dims_phase_complex);

#ifdef _USE_APERTURE_
    /* ----------------------------------------------
     * If aperture function available, get from root.
     * ----------------------------------------------
     */

        Array<precision> aperture(dims_phase_cropped);
        MPI_Bcast(aperture[0], aperture.get_size(), mpi_precision, 0, MPI_COMM_WORLD);
#endif

    /* -------------------------------
     * Import fft wisdom if available.
     * -------------------------------
     */

        fftw_import_wisdom_from_filename(io_t::rd_phs_wisdom_from.c_str());

    /*
     * Variable declaration:
     * --------------------------------
     * Name     Type        Description
     * --------------------------------
     * forward  fftw_plan   Re-usable FFTW plan for the forward transformation.
     */

        fftw_plan reverse = fftw_plan_dft_2d(dims_phase_complex[0], dims_phase_complex[1],\
                                             reinterpret_cast<fftw_complex*>(phase_fourier[0]),\
                                             reinterpret_cast<fftw_complex*>(phase_complex[0]),\
                                             FFTW_BACKWARD, FFTW_ESTIMATE);
   
    /*
     * Variable declaration:.
     * --------------------------------------------
     * Name                 Type        Description
     * --------------------------------------------
     * fried            precision   Fried parameter received from root.
     * piston           precision   Piston of the cropped phase-screen.
     */

        precision fried  = 0.0; 
        precision piston = 0.0;

    /* --------------------------------------------------
     * Enter loop to simulate phase-screens until killed.
     * --------------------------------------------------
     */

        MPI_Recv(&fried, 1,  mpi_precision, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
        while(mpi_status.MPI_TAG != mpi_cmds::kill){

            for(sizt ind = 0; ind < sims_t::realizations_per_fried; ind++){
            	    
            /* -------------------------------
             * Simulate a single phase-screen.
             * -------------------------------
             */

                make_phase_screen_fourier_shifted(phase_fourier, fried, sims_t::size_in_meters * aperture_t::sampling_factor);
                fftw_execute_dft(reverse, reinterpret_cast<fftw_complex*>(phase_fourier[0]), reinterpret_cast<fftw_complex*>(phase_complex[0]));

            /* --------------------------------------
             * Crop simulation to the requested size.
             * --------------------------------------
             */

                phase_cropped = phase_complex.get_crop(dims_crop_begin, dims_phase_cropped).cast_to_type<precision>();

#ifdef _USE_APERTURE_
            /* ---------------------------------------------------------------------------------
             * If aperture available, clip the phase-screen with the aperture and remove piston.
             * ---------------------------------------------------------------------------------
             */

                phase_cropped *= aperture;
                piston = phase_cropped.get_total() / aperture.get_total();
                phase_cropped -= aperture * piston;
#else
                piston = phase_cropped.get_total() / phase_cropped.get_size();
                phase_cropped -= piston;
#endif
                memcpy(phase_all[ind], phase_cropped[0], phase_cropped.get_size() * sizeof(precision));
            }
	    
        /* --------------------------------------------------------------
         * Send phase-screens to root, and receive a new fried parameter.
         * --------------------------------------------------------------
         */

            MPI_Send(phase_all[0], phase_all.get_size(), mpi_precision, 0, mpi_pmsg::ready, MPI_COMM_WORLD);
            MPI_Recv(&fried, 1,  mpi_precision, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
        }

    /* -------------------------
     * Write FFT wisdom to file.
     * -------------------------
     */
            
        fftw_export_wisdom_to_filename(io_t::rd_phs_wisdom_from.c_str());
        fftw_destroy_plan(reverse);
        fftw_cleanup();
    }

    MPI_Finalize();
    return(EXIT_SUCCESS);
}
