// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_CLASSIFIERS_FWD_ALG_CLASSIFIER_H
#define FLUOROSEQ_CLASSIFIERS_FWD_ALG_CLASSIFIER_H

#include <functional>

#include "classifiers/scored_classification.h"
#include "common/approximation_model.h"
#include "common/dye_seq.h"
#include "common/error_model.h"
#include "fwd_alg/binomial_transition.h"
#include "fwd_alg/detach_transition.h"
#include "fwd_alg/edman_transition.h"
#include "fwd_alg/emission.h"
#include "tensor/tensor.h"

namespace fluoroseq {

class FwdAlgClassifier {
public:
    FwdAlgClassifier(int num_timesteps,
                     int num_channels,
                     const ErrorModel& error_model,
                     const ApproximationModel& approximation_model,
                     int num_dye_seqs,
                     DyeSeq** dye_seqs,
                     int* dye_seqs_num_peptides,
                     int* dye_seqs_ids);
    ~FwdAlgClassifier();
    ScoredClassification classify(const Radiometry& radiometry);
    ScoredClassification* classify(int num_radiometries, 
                                    Radiometry** radiometries);

    DetachTransition* detach_transition;
    BinomialTransition* dud_transition;
    BinomialTransition* bleach_transition;
    std::function<double (double, int)> pdf;
    DyeSeq** dye_seqs;  // not owned
    int* dye_seqs_num_peptides;  // not owned
    int* dye_seqs_ids;  // not owned
    EdmanTransition** edman_transitions;
    Tensor** tensors;
    int num_dye_seqs;
    int num_timesteps;
    int num_channels;
    int max_num_dyes;
    int max_failed_edmans;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_CLASSIFIERS_FWD_ALG_CLASSIFIER_H
