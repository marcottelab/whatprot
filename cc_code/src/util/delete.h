/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_UTIL_DELETE_H
#define WHATPROT_UTIL_DELETE_H

namespace whatprot {

template <class T>
void delete_if_pointer(const T& t) {}

template <class T>
void delete_if_pointer(T* t) {
    delete t;
}

template <class T>
void delete_array(unsigned int size, T* array) {
    delete[] array;
}

template <class T>
void delete_array(unsigned int size, T** array) {
    for (unsigned int i = 0; i < size; i++) {
        delete array[i];
    }
    delete[] array;
}

}  // namespace whatprot

#endif  // WHATPROT_UTIL_DELETE_H
