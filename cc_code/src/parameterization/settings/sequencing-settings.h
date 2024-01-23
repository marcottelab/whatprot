/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_PARAMETERIZATION_SETTINGS_SEQUENCING_SETTINGS_H
#define WHATPROT_PARAMETERIZATION_SETTINGS_SEQUENCING_SETTINGS_H

namespace whatprot {

class SequencingSettings {
public:
    // We prune the emission matrix using this as a multiplier for the standard
    // deviation.
    double dist_cutoff;
    // We want to account for cross dye interactions when pruning the emission
    // matrix, but we need to set a reasonable limit, because some peptides
    // can potentially have enormous numbers of fluorophores. This maximum is
    // then applied to only the worst interaction (most FRET). Also self
    // interaction is handled much more easily.
    unsigned int cross_dye_maximum;
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_SETTINGS_SEQUENCING_SETTINGS_H