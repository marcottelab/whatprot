/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "hmm-classifier.h"

// Standard C++ library headers:
#include <functional>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/dye-track.h"
#include "common/scored-classification.h"
#include "hmm/hmm/peptide-hmm.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"
#include "util/range.h"

namespace whatprot {

namespace {
using std::function;
using std::vector;
}  // namespace

HMMClassifier::HMMClassifier(
        unsigned int num_timesteps,
        unsigned int num_channels,
        const SequencingModel& seq_model,
        const SequencingSettings& seq_settings,
        const vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs)
        : seq_model(seq_model),
          seq_settings(seq_settings),
          universal_precomputations(seq_model, num_channels),
          dye_seqs(dye_seqs),
          num_timesteps(num_timesteps),
          num_channels(num_channels) {
    max_num_dyes = 0;
    for (const SourcedData<DyeSeq, SourceCount<int>>& dye_seq : dye_seqs) {
        dye_seq_precomputations_vec.push_back(new DyeSeqPrecomputations(
                dye_seq.value, seq_model, num_timesteps, num_channels));
        const DyeSeqPrecomputations& back = *dye_seq_precomputations_vec.back();
        for (unsigned int c = 0; c < num_channels; c++) {
            // Two things to be aware of for the next line of code.
            //   * The first dimension of the tensor shape is always the
            //     timestep, so we need to add one to the channel to index to
            //     the correct dimension.
            //   * We subtract one from the tensor shape because the tensor
            //     shape for the channel is one larger than the number of dyes,
            //     as it needs to go from 0 to the number of dyes, inclusively.
            int num_dyes = back.tensor_shape[1 + c] - 1;
            if (num_dyes > max_num_dyes) {
                max_num_dyes = num_dyes;
            }
        }
    }
    universal_precomputations.set_max_num_dyes(max_num_dyes);
}

HMMClassifier::~HMMClassifier() {
    for (DyeSeqPrecomputations* ds_pre : dye_seq_precomputations_vec) {
        delete ds_pre;
    }
}

ScoredClassification HMMClassifier::classify(const Radiometry& radiometry) {
    return classify_helper<Range>(radiometry, Range(dye_seqs.size()));
}

ScoredClassification HMMClassifier::classify(
        const Radiometry& radiometry, const vector<int>& candidate_indices) {
    return classify_helper<const vector<int>&>(radiometry, candidate_indices);
}

vector<ScoredClassification> HMMClassifier::classify(
        const vector<Radiometry>& radiometries) {
    vector<ScoredClassification> results;
    results.resize(radiometries.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (unsigned int i = 0; i < radiometries.size(); i++) {
        results[i] = classify(radiometries[i]);
    }
    return results;
}

vector<double> HMMClassifier::score(const Radiometry& radiometry) {
    RadiometryPrecomputations radiometry_precomputations(
            radiometry, seq_model, seq_settings, max_num_dyes);
    vector<double> scores;
    for (unsigned int i = 0; i < dye_seq_precomputations_vec.size(); i++) {
        PeptideHMM hmm(num_timesteps,
                       num_channels,
                       *dye_seq_precomputations_vec[i],
                       radiometry_precomputations,
                       universal_precomputations);
        scores.push_back(hmm.probability());
    }
    return scores;
}

vector<vector<double>> HMMClassifier::score(
        const vector<Radiometry>& radiometries) {
    vector<vector<double>> all_scores;
    all_scores.resize(radiometries.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (unsigned int i = 0; i < radiometries.size(); i++) {
        all_scores[i] = score(radiometries[i]);
    }
    return all_scores;
}

}  // namespace whatprot
