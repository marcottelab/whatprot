/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_TENSOR_VECTOR_H
#define WHATPROT_TENSOR_VECTOR_H

namespace whatprot {

class Vector {
public:
    Vector(int length, int stride, double* values);
    double& operator[](int i);
    double operator[](int i) const;

    double* values;
    int length;
    int stride;
};

}  // namespace whatprot

#endif  // WHATPROT_TENSOR_VECTOR_H
