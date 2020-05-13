// Author: Matthew Beauregard Smith
#ifndef FLUOROSEQ_COMMON_APPROXIMATION_MODEL_H
#define FLUOROSEQ_COMMON_APPROXIMATION_MODEL_H

namespace fluoroseq {

class ApproximationModel {
public:
    ApproximationModel(int max_failed_edmans);

    int max_failed_edmans;
};

}  // fluoroseq

#endif  // FLUOROSEQ_COMMON_APPROXIMATION_MODEL_H