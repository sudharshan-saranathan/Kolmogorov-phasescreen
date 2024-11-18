#ifndef _LIBMEM_
#define _LIBMEM_

#include <cstdlib>

template <typename type>
class memory{
public:

/* 
 * Function declaration:
 * ----------------------------------------
 * Name         Return Type     Description
 * ----------------------------------------
 * allocate()   type*           Returns 1D pointer pointing to <src> with dimensions
 *                              xs. If <src> = nullptr, memory is allocated on the heap.
 */

    static type*    allocate(std::size_t xs, type* src = nullptr){
        type* data = src == nullptr ? new type[xs] : src;
        return(data);
    }

/* 
 * Function declaration:
 * ----------------------------------------
 * Name         Return Type     Description
 * ----------------------------------------
 * allocate()   type**          Returns a 2D pointer pointing to <src> with dimensions
 *                              (xs, ys). If <src> = nullptr, memory is allocated on the
 *                              heap.
 */

    static type**   allocate(std::size_t xs, std::size_t ys, type* src = nullptr){
        type** data = new type*[xs];
        if(data == nullptr)
            return nullptr;

        data[0] = src == nullptr ? new type[xs*ys]() : src;
        if(data[0] == nullptr){
            delete[] data;
            return nullptr;
        }

        for(std::size_t i=1; i<xs; i++)
            data[i] = data[i-1] + ys;
        
        return data;
    }

/* 
 * Function declaration:
 * ----------------------------------------
 * Name         Return Type     Description
 * ----------------------------------------
 * allocate()   type***         Returns a 3D pointer pointing to <src> with dimensions
 *                              (xs, ys, zs). If <src> = nullptr, memory is allocated on
 *                              the heap.
 */

    static type***  allocate(std::size_t xs, std::size_t ys, std::size_t zs, type* src = nullptr){
        type ***data = new type**[xs]();
        if(data == nullptr)
            return nullptr;

        data[0] = new type*[xs*ys]();
        if(data[0] == nullptr){
            delete[] data;
            return nullptr;
        }

        data[0][0] = src == nullptr ? new type[xs*ys*zs]() : src;
        if(data[0][0] == nullptr){
            delete[] data[0];
            delete[] data;
            return nullptr;
        }

        for(std::size_t i=0; i<xs; i++){
            *(data+i) = *data + i*ys;
            for(std::size_t j=0; j<ys; j++)
                *(*(data+i)+j) = **data + i*ys*zs + j*zs;
        }
        return data;
  }

/* 
 * Function declaration:
 * ----------------------------------------
 * Name         Return Type     Description
 * ----------------------------------------
 * allocate()   type****        Returns a 4D pointer pointing to <src> with dimensions
 *                              (xs, ys, zs, ws). If <src> = nullptr, memory is allocated
 *                              on the heap.
 */

    static type**** allocate(std::size_t xs, std::size_t ys, std::size_t zs, std::size_t ws, type* src = nullptr){
        type ****data = new type***[xs]();
        if(data==nullptr)
            return(nullptr);

        data[0] = new type**[xs*ys]();
        if(data[0]==nullptr){
            delete[] data;
            return nullptr;
        }

        data[0][0] = new type*[xs*ys*zs]();
        if(data[0][0]==nullptr){
            delete[] data[0];
            delete[] data;
            return(nullptr);
        }

        data[0][0][0] = src == nullptr ? new type[xs*ys*zs*ws]() : src;
        if(data[0][0][0] == nullptr){
            delete[] data[0][0];
            delete[] data[0];
            delete[] data;
            return nullptr;
        }

        for(std::size_t i=0; i<xs; i++){
            *(data+i) = *data + i*ys;
            for(std::size_t j=0; j<ys; j++){
                *(*(data+i)+j) = **data + i*ys*zs + j*zs;
                for(std::size_t k=0; k<zs; k++)
                    *(*(*(data+i)+j)+k) = ***data + i*ys*zs*ws + j*zs*ws + k*ws;
            }
        }
        return data;
    }

/* 
 * Function declaration:
 * ----------------------------------------
 * Name         Return Type     Description
 * ----------------------------------------
 * deallocate() type****        Frees a 4D pointer. If <free> = false, only
 *                              the higher dimensional pointers are freed.
 */

    static void deallocate(type**** data, bool owner = true){

        if(data != nullptr){
            if(data[0] != nullptr){
                if(data[0][0] != nullptr){
                    if(data[0][0][0] != nullptr){
                        if(owner)
                            delete[] data[0][0][0];
                        else
                            data[0][0][0] = nullptr;
                    }
                    delete[] data[0][0];
                }
                delete[] data[0];
            }
            delete[] data;
        }
    }

/* 
 * Function declaration:
 * ----------------------------------------
 * Name         Return Type     Description
 * ----------------------------------------
 * deallocate() type****        Frees a 3D pointer. If <free> = false, only
 *                              the higher dimensional pointers are freed.
 */

    static void deallocate(type***  data, bool owner = true){

        if(data != nullptr){
            if(data[0] != nullptr){
                if(data[0][0] != nullptr){
                    if(owner)
                        delete[] data[0][0];
                    else
                        data[0][0] = nullptr;
                }
                delete[] data[0];
            }
            delete[] data;
        }
    }

/* 
 * Function declaration:
 * ----------------------------------------
 * Name         Return Type     Description
 * ----------------------------------------
 * deallocate() type****        Frees a 2D pointer. If <free> = false, only 
 *                              the higher dimensional pointers are freed.
 */

    static void deallocate(type**   data, bool owner = true){
        if(data != nullptr){
            if(data[0] != nullptr){
                if(owner)
                    delete[] data[0];
                else
                    data[0] = nullptr;
            }
            delete[] data;
        }
    }

/* 
 * Function declaration:
 * ----------------------------------------
 * Name         Return Type     Description
 * ----------------------------------------
 * deallocate() type****        Frees a 1D pointer. If <free> = false, only 
 *                              the higher dimensional pointers are freed.
 */

    static void deallocate(type*    data, bool owner = true){
        if(data != nullptr){
            if(owner)
                delete[] data;
            else
                data = nullptr;
        }
    }
};

#endif
