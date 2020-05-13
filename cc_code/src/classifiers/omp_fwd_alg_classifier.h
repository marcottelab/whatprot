// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_CLASSIFIERS_OMP_FWD_ALG_CLASSIFIER_H
#define FLUOROSEQ_CLASSIFIERS_OMP_FWD_ALG_CLASSIFIER_H

#ifdef _OPENMP

#include "common/dye_seq.h"
#include "common/error_model.h"
#include "classifiers/fwd_alg_classifier.h"

namespace fluoroseq {

class OMPFwdAlgClassifier {
public:
    OMPFwdAlgClassifier(int num_timesteps,
                        int num_channels,
                        const ErrorModel& error_model,
                        const ApproximationModel& approximation_model,
                        int num_dye_seqs,
                        DyeSeq** dye_seqs);
    ~OMPFwdAlgClassifier();
    ScoredClassification classify(const Radiometry& radiometry);
    ScoredClassification* classify(int num_radiometries, 
                                    Radiometry** radiometries);

    FwdAlgClassifier** classifiers;
};

}  // namespace fluoroseq

#endif  // _OPENMP

#endif  // FLUOROSEQ_CLASSIFIERS_OMP_FWD_ALG_CLASSIFIER_H
