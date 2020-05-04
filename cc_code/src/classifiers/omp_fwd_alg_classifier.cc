// Author: Matthew Beauregard Smith (UT Austin)
#include "omp_fwd_alg_classifier.h"

#ifdef _OPENMP

#include <omp.h>

#include "classifiers/scored_classification.h"

namespace fluoroseq {

# pragma omp declare reduction(merge_scores                                    \
                               : ScoredClassification*                         \
                               : omp_out = merge_scores(omp_in, omp_out))      \
        initializer(omp_priv = new ScoredClassification())

OMPFwdAlgClassifier::OMPFwdAlgClassifier(int num_omp_dye_seq_groups,
                                         int num_timesteps,
                                         int num_channels,
                                         const ErrorModel& error_model,
                                         int num_dye_seqs,
                                         DyeSeq** dye_seqs)
        : num_omp_dye_seq_groups(num_omp_dye_seq_groups) {
    classifiers = new FwdAlgClassifier*[
            num_omp_dye_seq_groups * omp_get_max_threads()];
    #pragma omp parallel for
    for (int i = 0; i < omp_get_max_threads(); i++) {
        int start = 0;
        int end = num_dye_seqs / num_omp_dye_seq_groups;
        for (int j = 0; j < num_omp_dye_seq_groups; j++) {
            // Casting to long to avoid overflow. Probably unnecessary, but it
            // doesn't seem totally impossible that num_dye_seqs *
            // num_omp_dye_seq_groups could be greater than 4*10^9 (largest
            // value for an int).
            start = (long) num_dye_seqs * (long) j
                    / (long) num_omp_dye_seq_groups;
            end = (long) num_dye_seqs * (long) (j + 1)
                  / (long) num_omp_dye_seq_groups;
            classifiers[i * num_omp_dye_seq_groups + j] = new FwdAlgClassifier(
                    num_timesteps,
                    num_channels,
                    error_model,
                    end - start,
                    &dye_seqs[start]);
        }
    }
}

OMPFwdAlgClassifier::~OMPFwdAlgClassifier() {
    for (int i = 0; i < num_omp_dye_seq_groups * omp_get_max_threads(); i++) {
        delete classifiers[i];
    }
    delete[] classifiers;
}

ScoredClassification* OMPFwdAlgClassifier::classify(
        const Radiometry& radiometry) {
    ScoredClassification* best = new ScoredClassification();
    # pragma omp parallel for reduction(merge_scores : best)
    for (int i = 0; i < num_omp_dye_seq_groups; i++) {
        ScoredClassification* classification = classifiers[i]->classify(
                radiometry);
        best = merge_scores(best, classification);
    }
    return best;
}

}  // namespace flurooseq
#include <iostream>
using std::cout;
namespace fluoroseq {

ScoredClassification** OMPFwdAlgClassifier::classify(
        int num_radiometries, Radiometry** radiometries) {
    ScoredClassification** best = new ScoredClassification*[num_radiometries];
    for (int i = 0; i < num_radiometries; i++) {
        best[i] = new ScoredClassification();
    }
    omp_lock_t* lock_best = new omp_lock_t[num_radiometries];
    for (int i = 0; i < num_radiometries; i++) {
        omp_init_lock(&lock_best[i]);
    }
    #pragma omp parallel for
    for (int k = 0; k < num_omp_dye_seq_groups * num_radiometries; k++) {
        int i = k / num_radiometries;  // index of the omp dye seq group.
        int j = k % num_radiometries;  // index of the radiometry.
        ScoredClassification* classification = classifiers[
                omp_get_thread_num() * num_omp_dye_seq_groups + i]->classify(
                        *radiometries[j]);
        omp_set_lock(&lock_best[j]);
        best[j] = merge_scores(best[j], classification);
        omp_unset_lock(&lock_best[j]);
    }
    for (int i = 0; i < num_radiometries; i++) {
        omp_destroy_lock(&lock_best[i]);
    }
    delete[] lock_best;
    return best;
}

}  // namespace fluoroseq

#endif  // _OPENMP
