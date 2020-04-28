// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_COMMON_ERROR_MODEL_H
#define FLUOROSEQ_COMMON_ERROR_MODEL_H

namespace fluoroseq {

enum DistributionType {
// NORMAL,
LOGNORMAL,
};

class ErrorModel {
public:
    ErrorModel(double p_edman_failure,
               double p_detach,
               double p_bleach,
               double p_dud,
               DistributionType distribution_type,
               double mu,
               double sigma);
    double p_edman_failure;
    double p_detach;
    double p_bleach;
    double p_dud;
    DistributionType distribution_type;
    double mu;
    double sigma;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_COMMON_ERROR_MODEL_H
