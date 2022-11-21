/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_PARAMETERIZATION_MODEL_SEQUENCING_MODEL_H
#define WHATPROT_PARAMETERIZATION_MODEL_SEQUENCING_MODEL_H

// Standard C++ library headers:
#include <string>
#include <vector>

// Local project headers:
#include "parameterization/model/channel-model.h"

namespace whatprot {

class SequencingModel {
public:
    SequencingModel();
    SequencingModel(const std::string& seq_model_filename);
    SequencingModel(const SequencingModel& other);
    SequencingModel& operator=(const SequencingModel& other);
    SequencingModel(SequencingModel&& other);
    ~SequencingModel();
    SequencingModel with_mu_as_one() const;
    double relative_distance(const SequencingModel& sequencing_model) const;
    std::string debug_string() const;

    double p_edman_failure;
    double p_detach;
    std::vector<ChannelModel*> channel_models;
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_MODEL_SEQUENCING_MODEL_H
