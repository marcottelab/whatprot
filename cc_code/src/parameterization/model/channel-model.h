/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_PARAMETERIZATION_MODEL_CHANNEL_MODEL_H
#define WHATPROT_PARAMETERIZATION_MODEL_CHANNEL_MODEL_H

// Standard C++ library headers:
#include <cmath>
#include <functional>
#include <numbers>
#include <string>
#include <vector>

namespace whatprot {

class ChannelModel {
public:
    ChannelModel(unsigned int channel, unsigned int num_channels);
    // A virtual destructor is required since we have virtual methods, even
    // though we don't need special logic for anything.
    virtual ~ChannelModel();
    ChannelModel with_mu_as_one() const;
    // Note: This is virtual so that it can be mocked in tests.
    virtual double sigma(double amu) const;
    double distance(const ChannelModel& channel_model) const;
    std::string debug_string() const;

    // Note: N must be a numeric type, and the pointer must be indexable from 0
    // to num_channels - 1.
    template <class N>
    double adjusted_mu(const N* counts) const {
        double adjusted_mu = mu;
        adjusted_mu *= (double)counts[channel];
        for (unsigned int c = 0; c < num_channels; c++) {
            if (c == channel) {
                // For self quenching there is only an interaction if there are
                // two dyes on the channel.
                if (counts[c] > 1) {
                    adjusted_mu *= flat_interactions[c];
                    adjusted_mu *= std::pow(interactions[c], counts[c] - 1);
                }
            } else {
                if (counts[c] > 0) {
                    adjusted_mu *= flat_interactions[c];
                    adjusted_mu *= std::pow(interactions[c], counts[c]);
                }
            }
        }
        return adjusted_mu;
    }

    // Note: N must be a numeric type, and the pointer must be indexable from 0
    // to num_channels - 1.
    // Note: This needs to be virtual to allow mocking in tests, but template
    // functions can't be virtual. We work around this by making it a helper
    // function - so that we can still avoid code duplication.
    template <class N>
    double pdf_helper(double observed, const N* counts) const {
        static const double pi = 3.141592653589793238;
        double amu = adjusted_mu(counts);
        double offset = observed - amu;
        double s = sigma(amu);
        return (1.0 / (s * sqrt(2.0 * pi)))
               * exp(-offset * offset / (2.0 * s * s));
    }

    // virtual to allow mocking in tests.
    virtual double pdf(double observed, const unsigned int* counts) const;
    virtual double pdf(double observed, const short* counts) const;

    unsigned int channel;
    unsigned int num_channels;
    double p_bleach;
    double p_dud;
    double bg_sig;
    double mu;
    double sig;
    std::vector<double> interactions;
    std::vector<double> flat_interactions;
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_MODEL_CHANNEL_MODEL_H
