/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_TENSOR_VECTOR_H
#define FLUOROSEQ_TENSOR_VECTOR_H

namespace fluoroseq {

class Vector {
public:
    Vector(int length, int stride, double* values);
    double& operator[](int i);

    double* values;
    int length;
    int stride;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_TENSOR_VECTOR_H
