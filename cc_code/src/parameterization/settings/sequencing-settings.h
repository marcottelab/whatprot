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
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_SETTINGS_SEQUENCING_SETTINGS_H