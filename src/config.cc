#include "config.h"

#include <fstream>
#include <sstream>

string io_t::rd_image_from      = "image.fits";
string io_t::rd_fried_from      = "fried.fits";
string io_t::rd_basis_from      = "basis.fits";
string io_t::rd_coeff_from      = "coeff.fits";
string io_t::rd_aperture_from   = "pupil.fits";
string io_t::rd_psf_wisdom_from = "psf.wisdom";
string io_t::rd_phs_wisdom_from = "phs.wisdom";

string io_t::wr_psf_to      = "psf.fits";
string io_t::wr_phase_to    = "phase.fits";
string io_t::wr_image_to    = "image_convolved.fits";
string io_t::wr_residual_to = "residual.fits";

fmt_t format_t::rd_image    = fmt_t::FITS;
fmt_t format_t::rd_fried    = fmt_t::FITS;
fmt_t format_t::rd_basis    = fmt_t::FITS;
fmt_t format_t::rd_coeff    = fmt_t::FITS;
fmt_t format_t::rd_aperture = fmt_t::FITS;

fmt_t format_t::wr_psf      = fmt_t::FITS;
fmt_t format_t::wr_phase    = fmt_t::FITS;
fmt_t format_t::wr_image    = fmt_t::FITS;
fmt_t format_t::wr_residual = fmt_t::FITS;

bool   io_t::save    = true;
bool   io_t::clobber = false;

sizt  sims_t::realizations_per_fried = 400;
sizt  sims_t::size_x_in_pixels       = 92;
sizt  sims_t::size_y_in_pixels       = 92;
float sims_t::size_in_meters         = 10.0;

float aperture_t::size            = 1.0;
float aperture_t::sampling_factor = 1.5;
bool  aperture_t::make_airy_disk  = false;

float image_t::normalization     = 1600.0;
float image_t::original_sampling = 1.0;
float image_t::degraded_sampling = 1.0;

int config_parse(const char* filename){
  
    std::ifstream infile(filename);
    std::string line;
    if(!infile)
        return(EXIT_FAILURE);

    while(std::getline(infile, line)){
    
        std::stringstream tokens(line);
        std::string key, value;

        std::getline(tokens, key, '=');
        trim(key);

        std::getline(tokens, value, '=');
        tokens >> std::ws;
        trim(value);

        if(key == "image")
	        io_t::rd_image_from = value;

        else if(key == "fried")
	        io_t::rd_fried_from = value;

        else if(key == "basis")
	        io_t::rd_basis_from = value;

        else if(key == "aperture")
	        io_t::rd_aperture_from = value;

        else if(key == "coeff")
	        io_t::rd_coeff_from = value;
        
        else if(key == "fftw_psf")
	        io_t::rd_psf_wisdom_from = value;
        
        else if(key == "fftw_phase")
	        io_t::rd_phs_wisdom_from = value;
        
        else if(key == "phase")
	        io_t::wr_phase_to = value;
        
        else if(key == "convolved_images")
	        io_t::wr_image_to = value;
        
        else if(key == "residual")
	        io_t::wr_residual_to = value;
        
        else if(key == "psf"){
	        io_t::wr_psf_to = value;}
        
        else if(key == "realizations")
	        sims_t::realizations_per_fried = std::stoi(value);
        
        else if(key == "phase_size_x_in_pixels")
	        sims_t::size_x_in_pixels = std::stoi(value);
        
        else if(key == "phase_size_y_in_pixels")
	        sims_t::size_y_in_pixels = std::stoi(value);
        
        else if(key == "phase_size_in_meters")
	        sims_t::size_in_meters  = std::stof(value);
        
        else if(key == "aperture_size_in_meters")
	        aperture_t::size = std::stof(value);
        
        else if(key == "aperture_sampling")
	        aperture_t::sampling_factor = std::stof(value);
        
        else if(key == "airy_disk")
	        aperture_t::make_airy_disk = value == "true";
        
        else if(key == "normalization")
	        image_t::normalization = std::stof(value);

        else if(key == "original_sampling")
	        image_t::original_sampling = std::stof(value);
    
        else if(key == "degraded_sampling")
	        image_t::degraded_sampling = std::stof(value);
        
        else if(key == "save")
	        io_t::save = value == "true";
        
        else if(key == "clobber")
	        io_t::clobber = value == "true"; 
    
    }
    
    return(EXIT_SUCCESS);

}

int trim(std::string &str, char delimiter){

    sizt string_begin = str.find_first_not_of(delimiter);
    sizt string_end   = str.find_last_not_of (delimiter);

    if(string_begin != std::string::npos && string_end != std::string::npos)
        str = str.substr(string_begin, string_end - string_begin + 1);

    return(EXIT_SUCCESS);
}
