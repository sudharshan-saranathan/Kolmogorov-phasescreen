#include "fitsio.h"
#include "lib_mem.h"
#include "lib_array.h"

#include <ctime>
#include <cmath>
#include <string>
#include <cstring>
#include <fstream>
#include <unistd.h>
#include <exception>
#include <algorithm>

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   sizeof_vector(sizt_vector&, sizt index = 0)
 * Description:     Returns the product of elements in a vector
 * ------------------------------------------------------------
 */

sizt sizeof_vector(      sizt_vector &vector, sizt index){

#ifdef _CHKBNDS_
    if(index >= vector.size())
        throw std::logic_error("In function sizeof_vector(): Index out of bounds");
#endif

    sizt N = 1;
    for(sizt ind = index; ind < vector.size(); ind++){
        N *= vector[ind];
    }
    return(N);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   sizeof_vector(const sizt_vector&, sizt index = 0)
 * Description:     Returns the product of elements in a vector
 * ------------------------------------------------------------
 */

sizt sizeof_vector(const sizt_vector &vector, sizt index){

#ifdef _CHKBNDS_
    if(index >= vector.size())
        throw std::logic_error("In function sizeof_vector(): Index out of bounds");
#endif

    sizt N = 1;
    for(sizt ind = index; ind < vector.size(); ind++){
        N *= vector[ind];
    }
    return(N);
}

/* ------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::Array()
 * ----------------------------------------------
 */

template <class type>
Array<type>:: Array(){

/* ----------------------------------------------------
 * Nothing to do. Default values are already specified.
 * ---------------------------------------------------
 */

}

/* ------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::~Array()
 * ----------------------------------------------
 */

template <class type>
Array<type>::~Array(){

/* -----------------------
 * Free the data pointers.
 * -----------------------
 */

    memory<type>::deallocate(data_ptr_1D, this->owner);
    memory<type>::deallocate(data_ptr_2D, this->owner);
    memory<type>::deallocate(data_ptr_3D, this->owner);
    memory<type>::deallocate(data_ptr_4D, this->owner);

/* ------------------------
 * Reset flags, clear dims.
 * ------------------------
 */

    this->root_ptr    = nullptr;
    this->data_ptr_1D = nullptr;
    this->data_ptr_2D = nullptr;
    this->data_ptr_3D = nullptr;
    this->data_ptr_4D = nullptr;
    
    this->dims.clear();
    this->owner = false;
    this->stat  = false;
    this->size  = 0;
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::Array(const sizt_vector&, type*)
 * --------------------------------------------------------------
 */

template <class type>
Array<type>:: Array(const sizt_vector &dimensions, type* src){

/* ------------------------------
 * Copy <dimensions> into <dims>.
 * ------------------------------
 */

    this->dims  = dimensions;
    
/* ---------------------------------------------------------
 * Initialize <size> to the number of elements in the array.
 * ---------------------------------------------------------
 */
    
    this->size  = sizeof_vector(this->dims);

/* --------------------------------------------------------------------
 * If <src> is not a nullptr, the current instance acts like a pointer.
 * -------------------------------------------------------------------------
 */

    this->owner = src == nullptr ? true : false;

/* ---------------------------------------------------------
 * Initialize <data_ptr_...> according to <dims>, and <src>.
 * ---------------------------------------------------------
 */

    switch(this->dims.size()){
        case 1:  this->data_ptr_1D = memory<type>::allocate(this->dims[0], src);
                 this->stat        = this->data_ptr_1D == nullptr ? false : true;
                 this->root_ptr    = this->data_ptr_1D;
                 break;

        case 2:  this->data_ptr_2D = memory<type>::allocate(this->dims[0], this->dims[1], src);
                 this->stat        = this->data_ptr_2D == nullptr ? false : true;
                 this->root_ptr    = this->data_ptr_2D[0];
                 break;

        case 3:  this->data_ptr_3D = memory<type>::allocate(this->dims[0], this->dims[1], this->dims[2], src);
                 this->stat        = this->data_ptr_3D == nullptr ? false : true;
                 this->root_ptr    = this->data_ptr_3D[0][0];
                 break;

        case 4:  this->data_ptr_4D = memory<type>::allocate(this->dims[0], this->dims[1], this->dims[2], this->dims[3], src);
                 this->stat        = this->data_ptr_4D == nullptr ? false : true;
                 this->root_ptr    = this->data_ptr_4D[0][0][0];
                 break;

        default: throw std::logic_error("In function Array<type>::Array(), expected dimensions <= 4");
                 break;
    }
}

/* ------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::Array(const Array<type>&)
 * -------------------------------------------------------
 */

template <class type>
Array<type>:: Array(const Array<type> &src){

/* -----------------------------------------------------------------
 * Copy constructor always allocates new memory. Set <owner> = true.
 * -----------------------------------------------------------------
 */

    this->owner = true;

/* --------------------------------------------
 * Copy <src.nans>, <src.dims>, and <src.size>.
 * --------------------------------------------
 */

    this->nans  = src.nans;
    this->dims  = src.dims;
    this->size  = src.size;

/* ----------------------------------------------
 * Initialize <data_ptr_...> according to <dims>.
 * ----------------------------------------------
 */

    switch(this->dims.size()){
        case 1: this->data_ptr_1D = memory<type>::allocate(this->dims[0]);
                this->stat        = this->data_ptr_1D == nullptr ? false : true;
                this->root_ptr    = this->data_ptr_1D;
                break;

        case 2: this->data_ptr_2D = memory<type>::allocate(this->dims[0], this->dims[1]);
                this->stat        = this->data_ptr_2D == nullptr ? false : true;
                this->root_ptr    = this->data_ptr_2D[0];
                break;

        case 3: this->data_ptr_3D = memory<type>::allocate(this->dims[0], this->dims[1], this->dims[2]);
                this->stat        = this->data_ptr_3D == nullptr ? false : true;
                this->root_ptr    = this->data_ptr_3D[0][0];
                break;

        case 4: this->data_ptr_4D = memory<type>::allocate(this->dims[0], this->dims[1], this->dims[2], this->dims[3]);
                this->stat        = this->data_ptr_4D == nullptr ? false : true;
                this->root_ptr    = this->data_ptr_4D[0][0][0];
                break;
    }

/* ----------------------------------------------
 * Copy data from <src.root_ptr> into <root_ptr>.
 * ----------------------------------------------
 */

    if(this->stat == true && src.stat == true)
        std::memcpy(this->root_ptr, src.root_ptr, this->size * sizeof(type));

}

/* ------------------------------------
 * Function header:         lib_array.h
 * Function name:           Array<type>::get_owner()
 * -------------------------------------------------
 */

template <class type>
bool         Array<type>:: get_owner(){
    return(this->owner);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::get_stat()
 * ----------------------------------------
 */

template <class type>
bool         Array<type>:: get_stat(){
    return(this->stat);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::get_size()
 * ----------------------------------------
 */

template <class type>
sizt         Array<type>:: get_size(){
    return(this->size);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::get_dims(sizt)
 * --------------------------------------------
 */

template <class type>
sizt         Array<type>:: get_dims(sizt xs){

#ifdef _CHKBNDS_
    if(xs >= dims.size())
        throw std::runtime_error("In function Array<type>::get_dims(), index out of bounds");
#endif

        return(dims[xs]);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::get_stat()
 * ----------------------------------------
 */

template <class type>
sizt_vector  Array<type>:: get_dims(){
    return(this->dims);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::get_total()
 * -----------------------------------------
 */

template <class type>
type	     Array<type>:: get_total(){

#ifdef _CHKSTAT_
    if(this->stat == false)
        throw std::runtime_error("In function Array<type>::get_total(), cannot compute total of empty array");
#endif

    type total = static_cast<type>(0);
    for(sizt ind = 0; ind < this->size; ind++){
	    total += this->root_ptr[ind];
    }
    return(total);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator[](const sizt)
 * ----------------------------------------------------
 */

template <class type>
type* Array<type>::operator[](const sizt xs){

#ifdef _CHKSTAT_
    if(this->stat == false)
        return(nullptr);
#endif

#ifdef _CHKBNDS_
    if(xs >= this->dims[0])
       return(nullptr);
#endif

/* -------------------------------------
 * Return pointer to the slice at index.
 * -------------------------------------
 */

    switch(this->dims.size()){
        case 1:  return(this->data_ptr_1D + xs);
        case 2:  return(this->data_ptr_2D[xs]);
        case 3:  return(this->data_ptr_3D[xs][0]);
        case 4:  return(this->data_ptr_4D[xs][0][0]);
        default: return(nullptr);
    }
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator()(const sizt)
 * ----------------------------------------------------
 */

template <class type>
type& Array<type>::operator()(const sizt xs){
    
#ifdef _CHKSTAT_
    if(this->stat == false)
        throw std::runtime_error("In Array<type>::operator()(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims.size() != 1)
        throw std::runtime_error("In Array<type>::operator()(), expected " + std::to_string(this->dims.size()) + " argument(s)");
#endif

#ifdef _CHKBNDS_
    if(xs >= this->dims[0])
        throw std::range_error("In function Array<type>::operator()(), index out of bounds");
#endif

    return(this->data_ptr_1D[xs]);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator()(const sizt, const sizt)
 * ----------------------------------------------------------------
 */

template <class type>
type& Array<type>::operator()(const sizt xs, const sizt ys){

#ifdef _CHKSTAT_
    if(this->stat == false)
        throw std::runtime_error("In Array<type>::operator()(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims.size() != 2)
        throw std::runtime_error("In Array<type>::operator()(), expected " + std::to_string(this->dims.size()) + " argument(s)");
#endif

#ifdef _CHKBNDS_
    if(xs >= this->dims[0] || ys >= this->dims[1])
        throw std::range_error("In function Array<type>::operator()(), index out of bounds");
#endif

    return(this->data_ptr_2D[xs][ys]);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator()(const sizt, const sizt, const sizt)
 * ----------------------------------------------------------------------------
 */

template <class type>
type& Array<type>::operator()(const sizt xs, const sizt ys, const sizt zs){
    
#ifdef _CHKSTAT_
    if(this->stat == false)
        throw std::runtime_error("In Array<type>::operator()(), cannot perform operation on empty array");
#endif
    
#ifdef _CHKDIMS_
    if(this->dims.size() != 3)
        throw std::runtime_error("In Array<type>::operator()(), expected " + std::to_string(this->dims.size()) + " argument(s)");
#endif

#ifdef _CHKBNDS_
    if(xs >= this->dims[0] || ys >= this->dims[1] || zs >= this->dims[2])
        throw std::range_error("In function Array<type>::operator()(), index out of bounds");
#endif

    return(this->data_ptr_3D[xs][ys][zs]);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator()(const sizt, const sizt, const sizt, const sizt)
 * ----------------------------------------------------------------------------------------
 */

template <class type>
type& Array<type>::operator()(const sizt xs, const sizt ys, const sizt zs, const sizt ws){

#ifdef _CHKSTAT_
    if(this->stat == false)
        throw std::runtime_error("In Array<type>::operator()(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims.size() != 4)
        throw std::runtime_error("In Array<type>::operator()(), expected " + std::to_string(this->dims.size()) + " argument(s)");
#endif

#ifdef _CHKBNDS_
    if(xs >= this->dims[0] || ys >= this->dims[1] || zs >= this->dims[2] || ws >this->dims[3])
        throw std::range_error("In function Array<type>::operator()(), index out of bounds");
#endif

    return(this->data_ptr_4D[xs][ys][zs][ws]);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::swap(Array<type>&, Array<type>&)
 * --------------------------------------------------------------
 */

template <typename type>
void swap(Array<type>& src_a, Array<type>& src_b){

/* -----------------------
 * Swap all class methods.
 * -----------------------
 */

    std::swap(src_a.stat,  src_b.stat);
    std::swap(src_a.size,  src_b.size);
    std::swap(src_a.dims,  src_b.dims);
    std::swap(src_a.owner, src_b.owner);

    std::swap(src_a.root_ptr   , src_b.root_ptr);
    std::swap(src_a.data_ptr_1D, src_b.data_ptr_1D);
    std::swap(src_a.data_ptr_2D, src_b.data_ptr_2D);
    std::swap(src_a.data_ptr_3D, src_b.data_ptr_3D);
    std::swap(src_a.data_ptr_4D, src_b.data_ptr_4D);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator= (Array<type>)
 * -----------------------------------------------------
 */

template <class type>
Array<type>& Array<type>::operator= (Array<type> src){

/* ----------------
 * Copy-swap idiom.
 * ----------------
 */

    swap(*this, src);
    return(*this);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator+ (const Array<type>)
 * -----------------------------------------------------------
 */

template <class type>
Array<type>  Array<type>::operator+ (const Array<type> &src){

#ifdef _CHKSTAT_
    if(!this->stat || !src.stat)
        throw std::logic_error("In function Array<type>::operator+(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims != src.dims)
        throw std::logic_error("In function Array<type>::operator+(), expected argument with equal dimensions");
#endif

    Array<type> sum(src.dims);
    for(sizt ind = 0; ind < src.size; ind++){
        sum.root_ptr[ind] = this->root_ptr[ind] + src.root_ptr[ind];
    }

    return(sum);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator- (const Array<type>)
 * -----------------------------------------------------------
 */

template <class type>
Array<type>  Array<type>::operator- (const Array<type> &src){

#ifdef _CHKSTAT_
    if(!this->stat || !src.stat)
        throw std::logic_error("In function Array<type>::operator-(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims != src.dims)
        throw std::logic_error("In function Array<type>::operator-(), expected argument with equal dimensions");
#endif

    Array<type> diff(src.dims);
    for(sizt ind = 0; ind < src.size; ind++){
        diff.root_ptr[ind] = this->root_ptr[ind] - src.root_ptr[ind];
    }

    return(diff);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator* (const Array<type>)
 * -----------------------------------------------------------
 */

template <class type>
Array<type>  Array<type>::operator* (const Array<type> &src){

#ifdef _CHK_STAT
    if(!this->stat || !src.stat)
        throw std::logic_error("In Array<type>::operator*(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims != src.dims)
        throw std::logic_error("In Array<type>::operator*(), expected argument with equal dimensions");
#endif

    Array<type> product(src.dims);
    for(sizt ind = 0; ind < src.size; ind++){
        product.root_ptr[ind] = this->root_ptr[ind] * src.root_ptr[ind];
    }

    return(product);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator/ (const Array<type>)
 * -----------------------------------------------------------
 */

template <class type>
Array<type>  Array<type>::operator/ (const Array<type> &src){

#ifdef _CHKSTAT_
    if(!this->stat || !src.stat)
        throw std::logic_error("In function Array<type>::operator/(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims != src.dims)
        throw std::logic_error("In function Array<type>::operator/(), expected argument with equal dimensions");
#endif

    Array<type> div(src.dims);
    type null = static_cast<type>(0);

/* -----------------------------------------
 * Handle divide-by-zero with std::infinity.
 * -----------------------------------------
 */

    for(sizt ind = 0; ind < src.size; ind++){
        div.root_ptr[ind] = src.root_ptr[ind] == null ? std::numeric_limits<type>::infinity() : this->root_ptr[ind] / src.root_ptr[ind];
    }
    return(div);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator+=(const Array<type>&)
 * -----------------------------------------------------------
 */

template <class type>
void         Array<type>::operator+=(const Array<type> &src){

#ifdef _CHKSTAT_    
    if(!this->stat || !src.stat)
        throw std::logic_error("In function Array<type>::operator+=(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims != src.dims)
        throw std::logic_error("In function Array<type>::operator+=(), expected argument with equal dimensions");
#endif

    for(sizt ind = 0; ind < src.size; ind++){
        this->root_ptr[ind] += src.root_ptr[ind];
    }
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator-=(const Array<type>&)
 * -----------------------------------------------------------
 */

template <class type>
void         Array<type>::operator-=(const Array<type> &src){
    
#ifdef _CHKSTAT_    
    if(!this->stat || !src.stat)
        throw std::logic_error("In function Array<type>::operator-=(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims != src.dims)
        throw std::logic_error("In function Array<type>::operator-=(), expected argument with equal dimensions");
#endif

    for(sizt ind = 0; ind < src.size; ind++){
        this->root_ptr[ind] -= src.root_ptr[ind];
    }
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator*=(const Array<type>&)
 * -----------------------------------------------------------
 */

template <class type>
void         Array<type>::operator*=(const Array<type> &src){
    
#ifdef _CHKSTAT_    
    if(!this->stat || !src.stat)
        throw std::logic_error("In function Array<type>::operator*=(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims != src.dims)
        throw std::logic_error("In function Array<type>::operator*=(), expected argument with equal dimensions");
#endif

    for(sizt ind = 0; ind < src.size; ind++){
        this->root_ptr[ind] *= src.root_ptr[ind];
    }
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator/=(const Array<type>&)
 * -----------------------------------------------------------
 */

template <class type>
void         Array<type>::operator/=(const Array<type> &src){
    
#ifdef _CHKSTAT_    
    if(!this->stat || !src.stat)
        throw std::logic_error("In function Array<type>::operator/=(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims != src.dims)
        throw std::logic_error("In function Array<type>::operator/=(), expected argument with equal dimensions");
#endif

/* -------------------------------------------
 * Handle divide-by-zero with std::infinity().
 * -------------------------------------------
 */

    for(sizt ind = 0; ind < src.size; ind++)
        this->root_ptr[ind] /= src.root_ptr[ind] == static_cast<type>(0) ? std::numeric_limits<type>::infinity() : src.root_ptr[ind];
    
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator+(type)
 * ---------------------------------------------
 */

template <class type>
Array<type>  Array<type>::operator +(type value){

#ifdef _CHKSTAT_    
    if(!this->stat)
        throw std::logic_error("In function Array<type>::operator+(), cannot perform operation on empty array");
#endif

    Array<type> sum(this->dims);
    for(sizt ind = 0; ind < this->size; ind++){
        sum.root_ptr[ind] = this->root_ptr[ind] + value;
    }
    return(sum);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator-(type)
 * ---------------------------------------------
 */

template <class type>
Array<type>  Array<type>::operator -(type value){

#ifdef _CHKSTAT_
    if(!this->stat)
        throw std::logic_error("In function Array<type>::operator-(), cannot perform operation on empty array");
#endif

    Array<type> diff(this->dims);
    for(sizt ind = 0; ind < this->size; ind++){
        diff.root_ptr[ind] = this->root_ptr[ind] - value;
    }
    return(diff);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator*(type)
 * ---------------------------------------------
 */

template <class type>
Array<type>  Array<type>::operator *(type value){

#ifdef _CHKSTAT_
    if(!this->stat)
        throw std::logic_error("In function Array<type>::operator*(), cannot perform operation on empty array");
#endif

    Array<type> product(this->dims);
    for(sizt ind = 0; ind < this->size; ind++){
        product.root_ptr[ind] = this->root_ptr[ind] * value;
    }
    return(product);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator/(type)
 * ---------------------------------------------
 */

template <class type>
Array<type>  Array<type>::operator /(type value){

#ifdef _CHKSTAT_
    if(!this->stat)
        throw std::logic_error("In function Array<type>::operator/(), cannot perform operation on empty array");
#endif

    Array<type> div(this->dims);
    type null = static_cast<type>(0);

    if(value == null){
        for(sizt ind = 0; ind < this->size; ind++)
            div.root_ptr[ind] = std::numeric_limits<type>::infinity();

    }else{
        for(sizt ind = 0; ind < this->size; ind++)
            div.root_ptr[ind] = this->root_ptr[ind] / value;

    }
    
    return(div);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator+=(type)
 * ----------------------------------------------
 */

template <class type>
void         Array<type>::operator+=(type value){

#ifdef _CHKSTAT_
    if(!this->stat)
        throw std::logic_error("In function Array<type>::operator+=(), cannot perform operation on empty array");
#endif

    for(sizt ind = 0; ind < this->size; ind++){
        this->root_ptr[ind] += value;
    }
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator-=(type)
 * ----------------------------------------------
 */

template <class type>
void         Array<type>::operator-=(type value){

#ifdef _CHKSTAT_
    if(!this->stat)
        throw std::logic_error("In function Array<type>::operator-=(), cannot perform operation on empty array");
#endif

    for(sizt ind = 0; ind < this->size; ind++){
        this->root_ptr[ind] -= value;
    }
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator*=(type)
 * ---------------------------------------------
 */

template <class type>
void         Array<type>::operator*=(type value){

#ifdef _CHKSTAT_
    if(!this->stat)
        throw std::logic_error("In function Array<type>::operator*=(), cannot perform operation on empty array");
#endif

    for(sizt ind = 0; ind < this->size; ind++){
        this->root_ptr[ind] *= value;
    }
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::operator/=(type)
 * ----------------------------------------------
 */

template <class type>
void         Array<type>::operator/=(type value){

#ifdef _CHKSTAT_
    if(!this->stat)
        throw std::logic_error("In function Array<type>::operator/=(), cannot perform operation on empty array");
#endif

    type null = static_cast<type>(0);
    if(value == null){
        for(sizt ind = 0; ind < this->size; ind++)
            this->root_ptr[ind] = std::numeric_limits<type>::infinity();
    
    }else{
        for(sizt ind  = 0; ind < this->size; ind++)
            this->root_ptr[ind] /= value;

    }
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::get_slice(sizt, bool)
 * ---------------------------------------------------
 */

template <class type>
Array<type>  Array<type>::get_slice(sizt index, bool copy){

#ifdef _CHKSTAT_
    if(!this->stat)
        throw std::runtime_error("In function Array<type>::get_slice(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(index >= dims[0])
        throw std::range_error("In function Array<type>::get_slice(), index out of bounds");
    if(this->dims.size() == 1)
        throw std::range_error("In function Array<type>::get_slice(), cannot slice 1D array");
#endif

    type *slice_ptr = nullptr;
    switch(this->dims.size()){
        case 2: slice_ptr = this->data_ptr_2D[index];
                break;
        case 3: slice_ptr = this->data_ptr_3D[index][0];
                break;
        case 4: slice_ptr = this->data_ptr_4D[index][0][0];
                break;
        default:throw std::runtime_error("In function Array<type>::get_slice(), expected to slice a 2D/3D/4D array");
                break;
    }

    sizt_vector dims_slice(this->dims.begin() + 1, this->dims.end());
    Array<type> data_slice(dims_slice, copy == true ? nullptr : slice_ptr);
    
    if(copy == true)
        memcpy(data_slice[0], slice_ptr, sizeof_vector(dims_slice) * sizeof(type));

    return(data_slice);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::get_abs()
 * ---------------------------------------
 */

template <class type>
Array<type>  Array<type>::get_abs(){

#ifdef _CHKSTAT_
    if(this->stat == false)
        throw std::runtime_error("In function Array<type>::get_abs(), cannot find absolute value of empty array");
#endif

    Array<type> absolute(this->dims);
    for(sizt ind = 0; ind < this->size; ind++)
        absolute.root_ptr[ind] = std::abs(this->root_ptr[ind]);

    return(absolute);

}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::get_norm()
 * ----------------------------------------
 */

template <class type>
Array<type>  Array<type>::get_norm(){

#ifdef _CHKSTAT_
    if(this->stat == false)
        throw std::runtime_error("In function Array<type>::get_norm(), cannot find absolute value of empty array");
#endif

    Array<type> norm(this->dims);
    for(sizt ind = 0; ind < this->size; ind++)
        norm.root_ptr[ind] = std::abs(this->root_ptr[ind]) * std::abs(this->root_ptr[ind]);

    return(norm);

}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::get_pad(sizt_vector, sizt_vector, type)
 * ---------------------------------------------------------------------
 */

template <class type>
Array<type>  Array<type>::get_pad(sizt_vector dims_start, sizt_vector dims_pad, type pad_value){

#ifdef _CHKSTAT_
    if(this->stat == false)
        throw std::runtime_error("In function Array<type>::get_pad(), cannot pad empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims.size() != dims_start.size() || this->dims.size() != dims_pad.size())
        throw std::runtime_error("In function Array<type>::get_pad(), expected " + std::to_string(this->dims.size()) + "D vector(s) as argument(s)");
    for(sizt ind = 0; ind < this->dims.size(); ind++){
        if(this->dims[ind] + dims_start[ind] > dims_pad[ind])
            throw std::runtime_error("In function Array<type>::get_pad(), dimensions of the padded array are too small");
    }
#endif
    
    Array<type> array_padded(dims_pad);
    switch(this->dims.size()){

        case 1: for(sizt xpix = 0; xpix < this->dims[0]; xpix++)
                    array_padded.data_ptr_1D[xpix + dims_start[0]] = this->data_ptr_1D[xpix];

                break;

        case 2: for(sizt xpix = 0; xpix < this->dims[0]; xpix++)
                    for(sizt ypix = 0; ypix < this->dims[1]; ypix++)
                        array_padded(xpix + dims_start[0], ypix + dims_start[1]) = this->data_ptr_2D[xpix][ypix];
                
                break;

        case 3: for(sizt xpix = 0; xpix < this->dims[0]; xpix++)
                    for(sizt ypix = 0; ypix < this->dims[1]; ypix++)
                        for(sizt zpix = 0; zpix < this->dims[2]; zpix++)
                            array_padded(xpix + dims_start[0], ypix + dims_start[1], zpix + dims_start[2]) = this->data_ptr_3D[xpix][ypix][zpix];
                
                break;

        case 4: for(sizt xpix = 0; xpix < this->dims[0]; xpix++)
                    for(sizt ypix = 0; ypix < this->dims[1]; ypix++)
                        for(sizt zpix = 0; zpix < this->dims[2]; zpix++)
                            for(sizt wpix = 0; wpix < this->dims[3]; wpix++)
                                array_padded(xpix + dims_start[0], ypix + dims_start[1], zpix + dims_start[2], wpix + dims_start[3]) = this->data_ptr_4D[xpix][ypix][zpix][wpix];
                
                break;
    }
    
    return(array_padded);

}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::get_roll(sizt_vector, bool)
 * ---------------------------------------------------------
 */

template <class type>
Array<type>  Array<type>::get_roll(sizt_vector shift, bool clockwise){

#ifdef _CHKSTAT_
    if(this->stat == false)
        throw std::runtime_error("In function Array<type>::get_roll(), cannot find absolute value of empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims.size() != shift.size())
        throw std::runtime_error("In function Array<type>::get_roll(), empty argument");
#endif
    
    Array<type> array_rolled(this->dims);
    if(!clockwise){
        for(sizt ind = 0; ind < this->dims.size(); ind++)
            shift[ind] += (this->dims[ind] % 2);
    }

    switch(shift.size()){

        case 1: for(sizt xpix = 0; xpix < this->dims[0]; xpix++){
                    array_rolled.data_ptr_1D[xpix] = this->data_ptr_1D[(xpix + shift[0]) % this->dims[0]];
                }
                break;

        case 2: for(sizt xpix = 0; xpix < this->dims[0]; xpix++){
                    for(sizt ypix = 0; ypix < this->dims[1]; ypix++){
                        array_rolled.data_ptr_2D[xpix][ypix] = this->data_ptr_2D[(xpix + shift[0]) % this->dims[0]]\
                                                                                [(ypix + shift[1]) % this->dims[1]];
                    }
                }
                break;
        
        case 3: for(sizt xpix = 0; xpix < this->dims[0]; xpix++){
                    for(sizt ypix = 0; ypix < this->dims[1]; ypix++){
                        for(sizt zpix = 0; zpix < this->dims[2]; zpix++){
                            array_rolled.data_ptr_3D[xpix][ypix][zpix] = this->data_ptr_3D[(xpix + shift[0]) % this->dims[0]]\
                                                                                          [(ypix + shift[1]) % this->dims[1]]\
                                                                                          [(zpix + shift[2]) % this->dims[2]];
                        }
                    }
                }
                break;
        
        case 4: for(sizt xpix = 0; xpix < this->dims[0]; xpix++){
                    for(sizt ypix = 0; ypix < this->dims[1]; ypix++){
                        for(sizt zpix = 0; zpix < this->dims[2]; zpix++){
                            for(sizt wpix = 0; wpix < this->dims[3]; wpix++){
                                array_rolled.data_ptr_4D[xpix][ypix][zpix][wpix] = this->data_ptr_4D[(xpix + shift[0]) % this->dims[0]]\
                                                                                                    [(ypix + shift[1]) % this->dims[1]]\
                                                                                                    [(zpix + shift[2]) % this->dims[2]]\
                                                                                                    [(wpix + shift[3]) % this->dims[3]];
                            }
                        }
                    }
                }
                break;
    }

    return(array_rolled);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::get_crop(sizt_vector, sizt_vector, bool)
 * ----------------------------------------------------------------------
 */

template <typename type>
Array<type>  Array<type>::get_crop(sizt_vector dims_start, sizt_vector dims_type, bool vector_type){

    sizt_vector dims_end(this->dims.size());
    sizt_vector dims_sub(this->dims.size());

#ifdef _CHKSTAT_
    if(!this->stat)
        throw std::runtime_error("In function Array<type>::get_crop(), cannot perform operation on empty array");
#endif

#ifdef _CHKDIMS_
    if(this->dims.size() != dims_start.size() || this->dims.size() != dims_type.size())
        throw std::runtime_error("In function Array<type>::get_crop(), expected " + std::to_string(this->dims.size()) + "D vector(s)");
#endif

#ifdef _CHKBNDS_
    if(vector_type){
        for(sizt ind = 0; ind < this->dims.size(); ind++){ 
            dims_sub[ind] = dims_type[ind];
            dims_end[ind] = dims_type[ind] + dims_start[ind];
            if(dims_end[ind] >= this->dims[ind])
                throw std::range_error("In function Array<type>::get_crop(), array out of bounds\n");
        }
    }else{
        for(sizt ind = 0; ind < this->dims.size(); ind++){
            dims_end[ind] = dims_type[ind];
            if(dims_start[ind] >= dims_end[ind])
                throw std::runtime_error("In function Array<type>::get_crop(), expected dims_start[" + std::to_string(ind) + "] <= dims_end[" + std::to_string(ind) + "]");
            else
                dims_sub[ind] = dims_end[ind] - dims_start[ind];
        }
    }
#endif

    Array<type> array_cropped(dims_sub);
    switch(this->dims.size()){

        case 1: for(sizt xpix = 0; xpix < dims_sub[0]; xpix++){
                    array_cropped(xpix) = this->data_ptr_1D[xpix + dims_start[0]];
                }
                break;

        case 2: for(sizt xpix = 0; xpix < dims_sub[0]; xpix++){
                    for(sizt ypix = 0; ypix < dims_sub[1]; ypix++){
                        array_cropped(xpix, ypix) = this->data_ptr_2D[xpix + dims_start[0]][ypix + dims_start[1]];
                    }
                }
                break;

        case 3: for(sizt xpix = 0; xpix < dims_sub[0]; xpix++){
                    for(sizt ypix = 0; ypix < dims_sub[1]; ypix++){
                        for(sizt zpix = 0; zpix < dims_sub[2]; zpix++){
                            array_cropped(xpix, ypix, zpix) = this->data_ptr_3D[xpix + dims_start[0]][ypix + dims_start[1]][zpix + dims_start[2]];                        
                        }
                    }
                }
                break;

        case 4: for(sizt xpix = 0; xpix < dims_sub[0]; xpix++){
                    for(sizt ypix = 0; ypix < dims_sub[1]; ypix++){
                        for(sizt zpix = 0; zpix < dims_sub[2]; zpix++){
                            for(sizt wpix = 0; wpix < dims_sub[3]; wpix++){
                                array_cropped(xpix, ypix, zpix, wpix) = this->data_ptr_4D[xpix + dims_start[0]][ypix + dims_start[1]][zpix + dims_start[2]][wpix + dims_start[3]];
                            }
                        }
                    }
                }
                break;
    }

    return(array_cropped);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::rd_bin(const char*)
 * -------------------------------------------------
 */

template <class type>
int          Array<type>::rd_bin(const char *filename){

    std::ifstream file(filename, std::ios::binary | std::ios::in | std::ios::ate);
    sizt filesize = file.tellg();
    file.seekg(0, std::ios::beg);

    if(this->stat == false || this->size != (filesize / sizeof(type))){
        file.close();
        return(EXIT_FAILURE);
    }

    file.read(reinterpret_cast<char*>(this->root_ptr), filesize);
    file.close();

    return(EXIT_SUCCESS);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::wr_bin(const char*, bool)
 * -------------------------------------------------------
 */

template <class type>
int          Array<type>::wr_bin(const char *filename, bool clobber){

    if(!clobber && !access(filename, F_OK))
        return(EXIT_FAILURE);

    std::ofstream file(filename, std::ios::binary | std::ios::out);
    file.write(reinterpret_cast<char*>(this->root_ptr), this->size * sizeof(type));
    file.close();

    return(EXIT_SUCCESS);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::rd_fits(const char*)
 * --------------------------------------------------
 */

template <class type>
int          Array<type>::rd_fits(const char *filename){

/* ------------------------------------------------------------------------
 * Determine the datatype and bitpix corresponding to the current instance.
 * ------------------------------------------------------------------------
 */

    int this_bitpix   = 0;
    int this_datatype = 0;

    if(std::is_same<type, float>::value){
        this_bitpix   = -32;
        this_datatype = TFLOAT;
    }else if(std::is_same<type, double>::value){
        this_bitpix   = -64;
        this_datatype = TDOUBLE;
    }else if(std::is_same<type, std::complex<float>>::value){
        this_bitpix   = -32;
        this_datatype = TCOMPLEX;
    }else if(std::is_same<type, std::complex<double>>::value){
        this_bitpix   = -64;
        this_datatype = TDBLCOMPLEX;
    }else{
        return(EXIT_FAILURE);
    }

/* ---------------------------
 * Create temporary variables.
 * ---------------------------
 */

    fitsfile *file = nullptr;
    sizt count = 1;
    int n_axis = 0;
    int status = 0;
    int bitpix = 0;

    sizt_vector fpix;
    sizt_vector dims;

    fits_open_file(&file, filename, READONLY, &status);
    if(status != 0)
        return(status);
 
/* ---------------------------
 * Probe FITS file dimensions.
 * ---------------------------
 */

    fits_get_img_dim (file, &n_axis, &status); dims.resize(n_axis);
    fits_get_img_size(file,  n_axis, (long int*)dims.data(), &status);
    fits_get_img_type(file, &bitpix, &status);
    if(status != 0)
        return(status);

/* ---------------------------------------------------------------------
 * In CFITSIO convention, dimensions of the array are stored in reverse.
 * ---------------------------------------------------------------------
 */

    std::reverse(dims.begin(), dims.end());

/* --------------------------------------
 * Create <data> to store the FITS array.
 * --------------------------------------
 */

    Array<type> data(dims);
    if(data.get_stat() == false)
        return(EXIT_FAILURE);

    count = sizeof_vector(dims);
    fpix.resize(n_axis); std::fill(fpix.begin(), fpix.end(), 1);

/* ---------------------------------------------------------------------
 * Read data from file. <this_datatype> determined earlier is used here.
 * ---------------------------------------------------------------------
 */

    fits_read_pix(file, this_datatype, (long int*)fpix.data(), count, nullptr, data[0], nullptr, &status);
    if(status != 0)
        return(status);

    *this = data;
    
    fits_close_file(file, &status);
    return(EXIT_SUCCESS);
}

/* ----------------------------
 * Function header: lib_array.h
 * Function name:   Array<type>::wr_fits/(const char*, bool)
 * ---------------------------------------------------------
 */

template <class type>
int 	     Array<type>::wr_fits(const char *name, bool clobber){

/* ------------------------------------------------------------------------
 * Determine the datatype and bitpix corresponding to the current instance.
 * ------------------------------------------------------------------------
 */

    int this_bitpix   = 0;
    int this_datatype = 0;

    if(std::is_same<type, float>::value){
        this_bitpix   = -32;
        this_datatype = TFLOAT;
    }else if(std::is_same<type, double>::value){
        this_bitpix   = -64;
        this_datatype = TDOUBLE;
    }else if(std::is_same<type, std::complex<float>>::value){
        this_bitpix   = -32;
        this_datatype = TCOMPLEX;
    }else if(std::is_same<type, std::complex<double>>::value){
        this_bitpix   = -64;
        this_datatype = TDBLCOMPLEX;
    }else{
        return(EXIT_FAILURE);
    }

/*
 * Variable declaration.
 * ----------------------------------------
 * Name	        Type            Description
 * ----------------------------------------
 * status	    int             CFITSIO routines use this to set errors.
 * fileptr	    fitsfile        FITS file pointer.
 * filename	    std::string     FITS file name.
 * dimensions	std::vector     Dimensions of the array, reversed.
 */

    int         status  = 0;
    fitsfile   *fileptr = nullptr;
    std::string filename(name);
    sizt_vector dimensions = this->dims;

    filename = clobber == true ? "!" + filename : filename;
    std::reverse(std::begin(dimensions), std::end(dimensions));

/* -------------------------
 * Create FITS file pointer.
 * -------------------------
 */

    fits_create_file(&fileptr, filename.c_str(), &status);
    if(status != 0)
        return(status);

/* ------------------
 * Create FITS image.
 * ------------------
 */

    fits_create_img(fileptr, this_bitpix, dimensions.size(), (long int*)dimensions.data(), &status);
    if(status != 0)
        return(status);
    
/* -------------------
 * Write data to file.
 * -------------------
 */

    fits_write_img(fileptr, this_datatype, 1, this->size, this->root_ptr, &status);
    if(status != 0)
        return(status);
    
    fits_close_file(fileptr, &status);
    return(status);
}

template class Array<int>;
template class Array<cmpx>;
template class Array<float>;
template class Array<double>;
