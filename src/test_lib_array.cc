#define _SIZE_ 32

#include "lib_array.h"
#include "lib_phase.h"

#include <map>
#include <cmath>
#include <cstdio>
#include <cstdlib>

int main(int argc, char *argv[]){

    system("clear");
    fprintf(stdout, "------------------------\n");
    fprintf(stdout, "Test module for libarray\n");
    fprintf(stdout, "------------------------\n");

    sizt_vector dims1D{2140};
    sizt_vector dims2D{_SIZE_, _SIZE_};
    sizt_vector dims3D{_SIZE_, _SIZE_, _SIZE_};
    sizt_vector dims4D{_SIZE_, _SIZE_, _SIZE_, _SIZE_};

    Array<double> array1D(dims1D);
    if(array1D.get_stat()){
        fprintf(stdout, "(Info):\tInitialized 1D instance\n");
        fprintf(stdout, "(Info):\tTesting class methods:\n");

        /*
        Testing class methods
        */
        fprintf(stdout, "\tget_stat():\t[%s]\n",  array1D.get_stat() == true ? "True" : "False");
        fprintf(stdout, "\tget_size():\t[%lu]\n", array1D.get_size());
        fprintf(stdout, "\tget_dims():\t[%ld]\n", array1D.get_dims()[0]);
        fprintf(stdout, "\trd_bin()  :\t[%s]\n",  array1D.rd_bin("fried.binary") == EXIT_SUCCESS ? "OK" : "Failed");
        fprintf(stdout, "\twr_bin()  :\t[%s]\n",  array1D.wr_bin("fried.binary.copy") == EXIT_SUCCESS ? "OK" : "Failed");
        /*
        Testing operator()
        */
        for(sizt xs = 0; xs < _SIZE_; xs++)
            array1D(xs) = sqrt(xs);

        fprintf(stdout, "\toperator():\t[%0.1lf]\n", array1D(_SIZE_ - 1));

        /*
        Testing copy constructor
        */
        Array<double> array1D_copy(array1D);
        fprintf(stdout, "\tcopyctor():\t[%s]\n", array1D_copy.get_stat() == true ? "OK" : "FAILED");

        /*
        Testing operator =
        */
        array1D_copy = array1D;
        fprintf(stdout, "\toperator =:\t[%s]\n", array1D_copy.get_stat() == true ? "OK" : "FAILED");

        /*
        Testing operator +
        */
        Array<double> array1D_add(array1D + array1D_copy);
        fprintf(stdout, "\toperator +:\t[%0.1lf]\n", array1D_add(_SIZE_ - 1));

        /*
        Testing operator -
        */
        Array<double> array1D_diff(array1D - array1D_copy);
        fprintf(stdout, "\toperator -:\t[%0.1lf]\n", array1D_diff(_SIZE_ - 1));
    
        /*
        Testing operator *
        */
        Array<double> array1D_product(array1D * array1D_copy);
        fprintf(stdout, "\toperator *:\t[%0.1lf]\n", array1D_product(_SIZE_ - 1));
    
        /*
        Testing operator /
        */
        Array<double> array1D_div(array1D / array1D_copy);
        fprintf(stdout, "\toperator /:\t[%0.1lf]\n", array1D_div(_SIZE_ - 1));

        /*
        Testing operator +=
        */
        array1D_add += array1D;
        fprintf(stdout, "\toperator+=:\t[%0.1lf]\n", array1D_add(_SIZE_ - 1));

        /*
        Testing operator -=
        */
        array1D_diff -= array1D;
        fprintf(stdout, "\toperator-=:\t[%0.1lf]\n", array1D_diff(_SIZE_ - 1));
    
        /*
        Testing operator *=
        */
        array1D_product *= array1D;
        fprintf(stdout, "\toperator*=:\t[%0.1lf]\n", array1D_product(_SIZE_ - 1));
    
        /*
        Testing operator /=
        */
        array1D_div /= array1D;
        fprintf(stdout, "\toperator/=:\t[%0.1lf]\n\n", array1D_div(_SIZE_ - 1));

    }

    Array<double> array2D(dims2D);
    if(array2D.get_stat()){
        fprintf(stdout, "(Info):\tInitialized 2D instance\n");
        fprintf(stdout, "(Info):\tTesting class methods:\n");

        /*
        Testing class methods
        */
        fprintf(stdout, "\tget_stat():\t[%s]\n",  array2D.get_stat() == true ? "True" : "False");
        fprintf(stdout, "\tget_size():\t[%lu]\n", array2D.get_size());
        fprintf(stdout, "\tget_dims():\t[%ld, %ld]\n", array2D.get_dims()[0], array2D.get_dims()[1]);
        //fprintf(stdout, "\trd_bin()  :\t[%s]\n",  array2D.rd_bin("fried.image") == EXIT_SUCCESS ? "OK" : "Failed");
        //fprintf(stdout, "\twr_bin()  :\t[%s]\n",  array2D.wr_bin("fried.image.copy") == EXIT_SUCCESS ? "OK" : "Failed");

        /*
        Testing operator()
        */
        for(sizt xs = 0; xs < _SIZE_; xs++){
            for(sizt ys = 0; ys < _SIZE_; ys++){
                array2D(xs, ys) = sqrt(xs*ys);
            }
        }

        fprintf(stdout, "\toperator():\t[%0.1lf]\n", array2D(12, 13));
    
        /*
        Testing copy constructor
        */
        Array<double> array2D_copy(array2D);
        fprintf(stdout, "\tcopyctor():\t[%s]\n", array2D_copy.get_stat() == true ? "OK" : "FAILED");

        /*
        Testing operator =
        */
        array2D_copy = array2D;
        fprintf(stdout, "\toperator =:\t[%s]\n", array2D_copy.get_stat() == true ? "OK" : "FAILED");

        /*
        Testing operator +
        */
        Array<double> array2D_add(array2D + array2D_copy);
        fprintf(stdout, "\toperator +:\t[%0.1lf]\n", array2D_add(_SIZE_ - 1, _SIZE_ - 1));

        /*
        Testing operator -
        */
        Array<double> array2D_diff(array2D - array2D_copy);
        fprintf(stdout, "\toperator -:\t[%0.1lf]\n", array2D_diff(_SIZE_ - 1, _SIZE_ - 1));
    
        /*
        Testing operator *
        */
        Array<double> array2D_product(array2D * array2D_copy);
        fprintf(stdout, "\toperator *:\t[%0.1lf]\n", array2D_product(_SIZE_ - 1, _SIZE_ - 1));
    
        /*
        Testing operator /
        */
        Array<double> array2D_div(array2D / array2D_copy);
        fprintf(stdout, "\toperator /:\t[%0.1lf]\n", array2D_div(_SIZE_ - 1, _SIZE_ - 1));

        /*
        Testing operator +=
        */
        array2D_add += array2D;
        fprintf(stdout, "\toperator+=:\t[%0.1lf]\n", array2D_add(_SIZE_ - 1, _SIZE_ - 1));

        /*
        Testing operator -=
        */
        array2D_diff -= array2D;
        fprintf(stdout, "\toperator-=:\t[%0.1lf]\n", array2D_diff(_SIZE_ - 1, _SIZE_ - 1));
    
        /*
        Testing operator *=
        */
        array2D_product *= array2D;
        fprintf(stdout, "\toperator*=:\t[%0.1lf]\n", array2D_product(_SIZE_ - 1, _SIZE_ - 1));
    
        /*
        Testing operator /=
        */
        array2D_div /= array2D;
        fprintf(stdout, "\toperator/=:\t[%0.1lf]\n\n", array2D_div(_SIZE_ - 1, _SIZE_ - 1));

    }

    Array<double> array3D(dims3D);
    if(array3D.get_stat()){
        fprintf(stdout, "(Info):\tInitialized 3D instance\n");
        fprintf(stdout, "(Info):\tTesting class methods:\n");

        /*
        Testing class methods
        */
        fprintf(stdout, "\tget_stat():\t[%s]\n",  array3D.get_stat() == true ? "True" : "False");
        fprintf(stdout, "\tget_size():\t[%lu]\n", array3D.get_size());
        fprintf(stdout, "\tget_dims():\t[%ld, %ld, %ld]\n", array3D.get_dims()[0], array3D.get_dims()[1], array3D.get_dims()[2]);

        /*
        Testing operator()
        */
        for(sizt xs = 0; xs < _SIZE_; xs++){
            for(sizt ys = 0; ys < _SIZE_; ys++){
                for(sizt zs = 0; zs < _SIZE_; zs++){
                    array3D(xs, ys, zs) = sqrt(xs*ys*zs);
                }
            }
        }

        fprintf(stdout, "\toperator():\t[%0.1lf]\n", array3D(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));

        /*
        Testing copy constructor
        */
        Array<double> array3D_copy(array3D);
        fprintf(stdout, "\tcopyctor():\t[%s]\n", array3D_copy.get_stat() == true ? "OK" : "FAILED");

        /*
        Testing operator =
        */
        array3D_copy = array3D;
        fprintf(stdout, "\toperator =:\t[%s]\n", array3D_copy.get_stat() == true ? "OK" : "FAILED");

        /*
        Testing operator +
        */
        Array<double> array3D_add(array3D + array3D_copy);
        fprintf(stdout, "\toperator +:\t[%0.1lf]\n", array3D_add(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));

        /*
        Testing operator -
        */
        Array<double> array3D_diff(array3D - array3D_copy);
        fprintf(stdout, "\toperator -:\t[%0.1lf]\n", array3D_diff(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));
    
        /*
        Testing operator *
        */
        Array<double> array3D_product(array3D * array3D_copy);
        fprintf(stdout, "\toperator *:\t[%0.1lf]\n", array3D_product(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));
    
        /*
        Testing operator /
        */
        Array<double> array3D_div(array3D / array3D_copy);
        fprintf(stdout, "\toperator /:\t[%0.1lf]\n", array3D_div(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));

        /*
        Testing operator +=
        */
        array3D_add += array3D;
        fprintf(stdout, "\toperator+=:\t[%0.1lf]\n", array3D_add(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));

        /*
        Testing operator -=
        */
        array3D_diff -= array3D;
        fprintf(stdout, "\toperator-=:\t[%0.1lf]\n", array3D_diff(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));
    
        /*
        Testing operator *=
        */
        array3D_product *= array3D;
        fprintf(stdout, "\toperator*=:\t[%0.1lf]\n", array3D_product(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));
    
        /*
        Testing operator /=
        */
        array3D_div /= array3D;
        fprintf(stdout, "\toperator/=:\t[%0.1lf]\n\n", array3D_div(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));

    }

    Array<double> array4D(dims4D);
    if(array4D.get_stat()){
        fprintf(stdout, "(Info):\tInitialized 4D instance\n");
        fprintf(stdout, "(Info):\tTesting class methods:\n");

        /*
        Testing class methods
        */
        fprintf(stdout, "\tget_stat():\t[%s]\n",  array4D.get_stat() == true ? "True" : "False");
        fprintf(stdout, "\tget_size():\t[%lu]\n", array4D.get_size());
        fprintf(stdout, "\tget_dims():\t[%ld, %ld, %ld, %ld]\n", array4D.get_dims()[0], array4D.get_dims()[1], array4D.get_dims()[2], array4D.get_dims()[3]);

        /*
        Testing operator()
        */
        for(sizt xs = 0; xs < _SIZE_; xs++){
            for(sizt ys = 0; ys < _SIZE_; ys++){
                for(sizt zs = 0; zs < _SIZE_; zs++){
                    for(sizt ws = 0; ws < _SIZE_; ws++){
                        array4D(xs, ys, zs, ws) = sqrt(xs*ys*zs*ws);
                    }
                }
            }
        }

        fprintf(stdout, "\toperator():\t[%0.1lf]\n", array4D(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));

        /*
        Testing copy constructor
        */
        Array<double> array4D_copy(array4D);
        fprintf(stdout, "\tcopyctor():\t[%s]\n", array4D_copy.get_stat() == true ? "OK" : "FAILED");

        /*
        Testing operator =
        */
        array4D_copy = array4D;
        fprintf(stdout, "\toperator =:\t[%s]\n", array4D_copy.get_stat() == true ? "OK" : "FAILED");

        /*
        Testing operator +
        */
        Array<double> array4D_add(array4D + array4D_copy);
        fprintf(stdout, "\toperator +:\t[%0.1lf]\n", array4D_add(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));

        /*
        Testing operator -
        */
        Array<double> array4D_diff(array4D - array4D_copy);
        fprintf(stdout, "\toperator -:\t[%0.1lf]\n", array4D_diff(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));
    
        /*
        Testing operator *
        */
        Array<double> array4D_product(array4D * array4D_copy);
        fprintf(stdout, "\toperator *:\t[%0.1lf]\n", array4D_product(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));
    
        /*
        Testing operator /
        */
        Array<double> array4D_div(array4D / array4D_copy);
        fprintf(stdout, "\toperator /:\t[%0.1lf]\n", array4D_div(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));

        /*
        Testing operator +=
        */
        array4D_add += array4D;
        fprintf(stdout, "\toperator+=:\t[%0.1lf]\n", array4D_add(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));

        /*
        Testing operator -=
        */
        array4D_diff -= array4D;
        fprintf(stdout, "\toperator-=:\t[%0.1lf]\n", array4D_diff(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));
    
        /*
        Testing operator *=
        */
        array4D_product *= array4D;
        fprintf(stdout, "\toperator*=:\t[%0.1lf]\n", array4D_product(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));
    
        /*
        Testing operator /=
        */
        array4D_div /= array4D;
        fprintf(stdout, "\toperator/=:\t[%0.1lf]\n\n", array4D_div(_SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1, _SIZE_ - 1));

    }

    return(EXIT_SUCCESS);
}
