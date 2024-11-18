#ifndef _LIBARRAY_
#define _LIBARRAY_

#define _USE_APERTURE_
#define _USE_SINGLE_PRECISION_

#define _CHKSTAT_
#define _CHKDIMS_
#define _CHKBNDS_

#include "config.h"
#include "fitsio.h"

#include <vector>
#include <cstdlib>
#include <typeinfo>
#include <stdexcept>

/*
 * Function declaration
 * ----------------------------------------------------------------
 * Name                                 Return type     Description 
 * ----------------------------------------------------------------
 * sizeof_vector(const sizt_vector&)    sizt            Returns the product of the vector elements, this
 *                                                      definition is used for vectors whose scope is
 *                                                      temporary.
 *
 * sizeof_vector(sizt_vector&)          sizt            Returns the product of the vector elements, this
 *                                                      definition is used when a valid object exists.
 */

sizt sizeof_vector(      sizt_vector&, sizt index = 0);
sizt sizeof_vector(const sizt_vector&, sizt index = 0);

/* --------------------
 * Template class Array
 * --------------------
 */

template <class type>
class Array{
private:

/*
 * Variable declaration (PRIVATE)
 * ------------------------------------
 * Name     Type            Description
 * ------------------------------------
 * dims     sizt_vector     A vector storing the dimensions of the array. The number of elements
 *                          in <dims> indicate the rank of the array (1D, 2D, 3D or 4D), while
 *                          each element in <dims> represents the size of the array along that
 *                          axis. For example, a 512 x 512 array would have <dims> = {512, 512}.
 *                      
 * owner    bool            True:   If <root_ptr> points to memory that was allocated on the heap
 *                                  by the current instance i.e. the memory block "belongs" to the
 *                                  current instance. The destructor frees the access pointers as
 *                                  well as the <root_ptr>.
 *                                  
 *                          False:  If <root_ptr> points to memory that was allocated on the heap
 *                                  by a different instance i.e. the current instance is like a 
 *                                  pointer. The destructor frees the access pointers, but not the
 *                                  <root_ptr> so as to prevent multiple frees of the same memory.
 *
 * nans     bool            True:   If the array has std::nans.
 *
 *                          False:  If the array does not have std::nans.
 *
 * stat     bool            True:   If <root_ptr> is pointing to a valid memory address. Always
 *                                  check if <stat> == true before accessing data to avoid seg-
 *                                  faults.
 *
 *                          False:  If <root_ptr> = nullptr. For example, when a new instance is
 *                                  initialized with the default constructor. 
 *
 * size     sizt            Number of elements in the array, obtained as the product of elements
 *                          in <dims>.
 */

    sizt_vector dims;
    bool owner = true;
    bool nans  = false;
    bool stat  = false;
    sizt size  = 1;

/*
 * Pointers declaration (PRIVATE)
 * ------------------------------------
 * Name         Type        Description
 * ------------------------------------
 * data_ptr_1D  type*       Pointer for accessing elements of a 1D array using matrix syntax.
 *                          For example, <data_ptr_1D[5]> accesses the 6th element.
 *
 * data_ptr_2D  type**      Pointers for accessing elements of a 2D array using matrix syntax.
 *                          For example, <data_ptr_2D[i][j]> accesses the element at row i, column j.
 *                          In pointer arithmetic, <data_ptr_2D[i][j]> = <root_ptr[i*dims[1] + j]>;
 *            
 * data_ptr_3D  type***     Pointers for accessing elements of a 3D array using matrix syntax.
 *                          See the description of <data_ptr_2D>.
 *
 * data_ptr_4D  type****    Pointers for accessing elements of a 4D array using matrix syntax.
 *                          See the description of <data_ptr_2D>.
 *
 * root_ptr     type*       Pointer pointing to the contiguous block of memory that stores all data.
 */

    type*    data_ptr_1D = nullptr;
    type**   data_ptr_2D = nullptr;
    type***  data_ptr_3D = nullptr;
    type**** data_ptr_4D = nullptr;
    type*    root_ptr    = nullptr;

public:

/*
 * Constructors / destructor declaration (PUBLIC)
 * ----------------------------------------------------
 * Name         Type                        Description
 * ----------------------------------------------------
 * Array()     ~Array()                     Destructor, uses <copy> to determine if <root_ptr>
 *                                          must also be freed.
 *
 * Array()      Array()                     Default constructor.
 *
 * Array()      Array(const sizt_vector&,\  Constructor that takes <dims> as the first argument.
 *                    type* src = nullptr)  If the second argument is nullptr, memory is allocated
 *                                          on the heap and stored in <root_ptr>. Otherwise, <src>
 *                                          is copied to <root_ptr>. 
 *
 *                                          Note: It is the user's responsibility to ensure that <src>
 *                                                can hold as many elements as sizeof_vector(<dims>).
 *
 * Array()      Array(const Array<type>&)   Copy constructor. Memory is newly allocated.
 */

   ~Array();
    Array();
    Array(const sizt_vector&, type* src = nullptr);
    Array(const Array<type>&);

public:

/* 
 * Class methods declaration (PUBLIC)
 * See lib_array.cc for definitions
 *
 * --------------------------------------------
 * Name             Return type     Description
 * --------------------------------------------
 * get_owner()      bool            Returns <owner>.
 * get_stat()       bool            Returns <stat>.
 * get_size()       sizt            Returns <size>.
 * get_total()      type            Returns the sum of all elements in the array.
 * get_dims(sizt)   sizt            Returns the size of the requested axis.
 * get_dims()       sizt_vector     Returns <dims>.
 */

    bool        get_owner();
    bool        get_stat ();
    sizt        get_size ();
    type	    get_total();
    sizt        get_dims (sizt);
    sizt_vector get_dims ();

public:

/*
 * Overloaded operators declaration (PUBLIC) 
 * See lib_array.cc for definitions
 * --------------------------------------------
 * Name             Return type     Description
 * --------------------------------------------
 * operator[]       type*           Returns a pointer to the memory location of the requested
 *                                  index in the first axis i.e. if the array is a 512 x 512
 *                                  image, array[100] returns a pointer to the starting of
 *                                  the 100th row. This function is useful for reading or
 *                                  writing slices of data at a time.
 * 
 * operator()       type&           Provides read / write access to the elements of the array,
 *                                  overloaded for 1D, 2D, 3D, and 4D arrays, throws exception
 *                                  if the index is out of bounds.
 *
 * operator=()      Array<type>&    Assignment operator for the class that used the copy-swap 
 *                                  idiom. 
 *
 * operator+()      Array<type>     Adds the argument to the array elements and returns the sum.
 *                                  If argument is an array, addition is performed element-wise
 *                                  i.e. dimensions must be the same.
 * 
 * operator-()      Array<type>     Subtracts the argument from the array elements and returns
 *                                  the difference. If argument is an array, subtraction is
 *                                  performed element-wise i.e. dimensions must be the
 *                                  same.
 *                                  
 * operator*()      Array<type>     Multiplies the array elements with the argument and returns
 *                                  the product. If argument is an array, multiplication is 
 *                                  performed element-wise i.e. dimensions must be the same. 
 *
 * operator/()      Array<type>     Divides the array elements with the argument and returns the
 *                                  quotient. If argument is an array, division is performed
 *                                  element-wise i.e. dimensions must be the same.
 *
 *                                  Note: If the divisor has zeros, the corresponding elements
 *                                        in the return array are set to std::nans.
 * 
 * operator+=()     void            Adds the argument to the array elements. If argument is an 
 *                                  array, addition is performed element-wise i.e. dimensions must
 *                                  match. 
 *                       
 * operator-=()     void            Subtracts the argument from the array elements. If argument
 *                                  is an array, subtraction is performed element-wise i.e. dimensions
 *                                  must match. 
 *
 * operator*=()     void            Multiplies the array elements with the argument. If argument
 *                                  is an array, multiplication is performed element-wise i.e.
 *                                  dimensions must match. 
 * 
 * operator/=()     void            Divides the array elements with the argument. If argument is
 *                                  an array, division is performed element-wise i.e. dimensions
 *                                  must match.
 *
 *                                  Note: If the divisor has zeros, the corresponding elements in 
 *                                        the array are set to std::nans.
 *
 * -----------------------------------
 * Usage tip for overloaded operators:
 * -----------------------------------
 *
 * Use the shorthand operators (+=, -=, *=, /=) if a new array is not required i.e. instead of:
 *  
 * array_1 = array_1 + array_2, (or)
 * array_1 = array_1 * 3
 *
 * use:
 * 
 * array_1 += array_2, (or)
 * array_1 *= 3
 *
 * This GINORMOUSLY saves memory usage (and time).
 */
  
    type* operator[](const sizt);
    type& operator()(const sizt);
    type& operator()(const sizt, const sizt);
    type& operator()(const sizt, const sizt, const sizt);
    type& operator()(const sizt, const sizt, const sizt, const sizt);
    
    Array<type>& operator = (      Array<type> );
    Array<type>  operator + (const Array<type>&);
    Array<type>  operator - (const Array<type>&);
    Array<type>  operator * (const Array<type>&);
    Array<type>  operator / (const Array<type>&);
    void         operator+= (const Array<type>&);
    void         operator-= (const Array<type>&);
    void         operator*= (const Array<type>&);
    void         operator/= (const Array<type>&);

    Array<type>  operator + (type);
    Array<type>  operator - (type);
    Array<type>  operator * (type);
    Array<type>  operator / (type);
    void         operator+= (type);
    void         operator-= (type);
    void         operator*= (type);
    void         operator/= (type);
 
public:

/* 
 * Class methods declaration (PUBLIC)
 * See lib_array.cc for definitions
 * --------------------------------------------
 * Name             Type            Description
 * --------------------------------------------
 * get_abs  ()      Array<type>     Returns the std::abs() of the array.
 * get_norm ()      Array<type>     Returns the std::norm() of the array.
 * get_pad  ()      Array<type>     Returns the array padded with zeros.
 * get_roll ()      Array<type>     Returns a cyclically rolled array.
 * get_crop ()      Array<type>     Returns a subset of the array.
 * get_slice()      Array<type>     Returns a slice of the array at the specified index.
 */

    Array<type> get_abs  ();
    Array<type> get_norm ();
    Array<type> get_pad  (sizt_vector, sizt_vector, type pad_value = static_cast<type>(0));
    Array<type> get_roll (sizt_vector, bool clockwise = true);
    Array<type> get_crop (sizt_vector, sizt_vector, bool vector_type = true);
    Array<type> get_slice(sizt, bool copy = true);

public:

/*
 * Class methods declaration (PUBLIC)
 * See lib_array.cc for definitions
 * --------------------------------------------
 * Name             Return type     Description
 * --------------------------------------------
 * swap()           void            Friend function to swap two instances. Required for the copy-swap idiom.
 *
 * rd_bin()         int             Read array from a binary file, the fastest option.
 *
 * wr_bin()         int             Write array to a binary file, the fastest option.
 *
 * cast_to_type()   void            Type-casting function for the class. Casts the elements
 *                                  of the array into the type of the argument, stores the type-casted
 *                                  values in the argument. Dimensions must match, throws exception if 
 *                                  they don't. See inline for definition.
 * 
 * cast_to_type()   Array<TYPE>     Type-casting function for the class. Returns a new type-casted 
 *                                  array.
 * 
 */

    template <typename TYPE>
    friend void swap(Array<TYPE>&, Array<TYPE>&);

    int rd_bin(const char*);
    int wr_bin(const char*,  bool clobber=false);

    template <typename TYPE>
    void cast_to_type(Array<TYPE>& array_to_cast_to){

    /* ---------------------------------------------
     * Check if dimensions of the input array match.
     * ---------------------------------------------
     */

        sizt_vector store_dims = array_to_cast_to.get_dims();
        if(this->dims.size() != store_dims.size())
            throw std::runtime_error("In function Array<type>::cast_to_type(), expected " + std::to_string(this->dims.size()) + "D argument");

        for(sizt ind = 0; ind < this->dims.size(); ind++){
            if(this->dims[ind] != store_dims[ind])
                throw std::runtime_error("In function Array<type>::cast_to_type(), expected dims[" + std::to_string(ind) + "] = " + std::to_string(this->dims[ind]) + " of argument");
        }

        if(this->stat == false)
            throw std::runtime_error("In function Array<type>::cast_to_type(), cannot cast an empty array");

    /* -----------------------
     * Cast the array to TYPE.
     * -----------------------
     * If the array being casted is of type std::complex, only the real part is casted.
     */ 

        for(sizt ind = 0; ind < this->size; ind++)
            *(array_to_cast_to[0] + ind) = static_cast<TYPE>(std::real(this->root_ptr[ind]));

    }

    template <typename TYPE>
    Array<TYPE> cast_to_type(){        

        Array<TYPE> array_to_cast_to(this->dims);
        if(this->stat == false)
            return array_to_cast_to;

    /* -----------------------
     * Cast the array to TYPE.
     * -----------------------
     * If the array being casted is of type std::complex<>, pick only the real part.
     */

        for(sizt ind = 0; ind < this->size; ind++)
            *(array_to_cast_to[0] + ind) = static_cast<TYPE>(std::real(this->root_ptr[ind]));
        
        return(array_to_cast_to);
    }

/*
 * Class methods declaration (PUBLIC)
 * See lib_array.cc for definitions
 * ------------------------------------
 * Name             Type    Description
 * ------------------------------------
 * rd_fits()        int     Read into array elements from FITS file, requires
 *                          the "cfitsio" library. Can be called on allocated
 *                          or empty instance. If read failed, fitsio error
 *                          codes are returned:
 *
 *                          https://heasarc.gsfc.nasa.gov/docs/software/fitsio/quick/node26.html
 *
 *                          A few frequently encountered error codes are listed below:
 *                          ----------------------------------------------------------
 *                          Err code = 101 (Input and output are the same file)
 *                          Err code = 104 (If the file could not be read. For example, it doesn't exist)
 *                          Err code = 108 (Error reading from FITS file)
 *                          ----------------------------------------------------------
 *
 * wr_fits()        int     Writes array elements to FITS file, requires the 
 *                          "cfitsio" library. If write failed, fitsio error
 *                          codes are returned:
 *
 *                          https://heasarc.gsfc.nasa.gov/docs/software/fitsio/quick/node26.html
 *
 *                          A few frequently encountered error codes are listed below:
 *                          ----------------------------------------------------------
 *                          Err code = 101 (Input and output are the same file)
 *                          Err code = 105 (If th file could not be created. For example, the filename includes a directory path that does not exist)
 *                          Err code = 106 (Error writing to FITS file)
 *                          ----------------------------------------------------------
 */

    int rd_fits(const char*);
    int wr_fits(const char*, bool clobber=false);
};

#endif
