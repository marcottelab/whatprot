// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_UTIL_DELETE_H
#define FLUOROSEQ_UTIL_DELETE_H

namespace fluoroseq {

template<class T>
void delete_if_pointer(const T& t) {}

template<class T>
void delete_if_pointer(T* t) {
    delete t;
}

template<class T> 
void delete_array(int size, T* array) {
    delete[] array;
}

template<class T>
void delete_array(int size, T** array) {
    for (int i = 0; i < size; i++) {
        delete array[i];
    }
    delete[] array;
}

}  // namespace fluoroseq

#endif  // FLUOROSEQ_UTIL_DELETE_H
