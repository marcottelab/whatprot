// Author: Matthew Beauregard Smith (UT Austin)
#include "error_model.h"

namespace fluoroseq {

ErrorModel::ErrorModel(double p_edman_failure,
                       double p_detach,
                       double p_bleach,
                       double p_dud,
                       DistributionType distribution_type,
                       double mu,
                       double sigma) : p_edman_failure(p_edman_failure),
                                       p_detach(p_detach),
                                       p_bleach(p_bleach),
                                       p_dud(p_dud),
                                       distribution_type(distribution_type),
                                       mu(mu),
                                       sigma(sigma) {}

}  // namespace fluoroseq
