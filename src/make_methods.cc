#include "lib_array.h"
#include "fitsio.h"
#include <cstdlib>
#include <cmath>

int main(){

    sizt_vector dims{64, 64, 64};
    sizt_vector dims_start{16, 16};
    sizt_vector dims_pad{128, 128};

    Array<double> array(dims);
    Array<double> array_padded(dims_pad);

    for(sizt xpix = 0; xpix < dims[0]; xpix++){
        for(sizt ypix = 0; ypix < dims[1]; ypix++){
            for(sizt zpix = 0; zpix < dims[2]; zpix++){
                array(xpix, ypix, zpix) = sqrt(xpix + ypix*zpix);
            }
        }
    }

    array_padded = array.get_slice(32).get_pad(dims_start, dims_pad);
    array_padded.wr_fits(std::string("array_padded.fits").c_str());

    return(0);
}
