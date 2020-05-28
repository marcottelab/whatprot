// Author: Matthew Beauregard Smith (UT Austin)
#include "omp_fwd_alg_classifier.h"

#ifdef _OPENMP

#include <omp.h>

#include "common/scored_classification.h"

namespace fluoroseq {

OMPFwdAlgClassifier::OMPFwdAlgClassifier(
        int num_timesteps,
        int num_channels,
        const ErrorModel& error_model,
        const ApproximationModel& approximation_model,
        int num_dye_seqs,
        DyeSeq** dye_seqs,
        int* dye_seqs_num_peptides,
        int* dye_seqs_ids) {
    classifiers = new FwdAlgClassifier*[omp_get_max_threads()];
    for (int i = 0; i < omp_get_max_threads(); i++) {
        classifiers[i] = new FwdAlgClassifier(num_timesteps,
                                              num_channels,
                                              error_model,
                                              approximation_model,
                                              num_dye_seqs,
                                              dye_seqs,
                                              dye_seqs_num_peptides,
                                              dye_seqs_ids);
    }
}

OMPFwdAlgClassifier::~OMPFwdAlgClassifier() {
    for (int i = 0; i < omp_get_max_threads(); i++) {
        delete classifiers[i];
    }
    delete[] classifiers;
}

ScoredClassification OMPFwdAlgClassifier::classify(
        const Radiometry& radiometry) {
    return classifiers[0]->classify(radiometry);
}

ScoredClassification* OMPFwdAlgClassifier::classify(
        int num_radiometries, Radiometry** radiometries) {
    ScoredClassification* result = new ScoredClassification[num_radiometries]();
    #pragma omp parallel for schedule(guided)
    for (int j = 0; j < num_radiometries; j++) {
        result[j] = classifiers[omp_get_thread_num()]->classify(
                *radiometries[j]);
    }
    return result;
}

}  // namespace fluoroseq

#endif  // _OPENMP
