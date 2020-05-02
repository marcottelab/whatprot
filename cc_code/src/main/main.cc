// Author: Matthew Beauregard Smith (UT Austin)
// Simple application to read in a TSV file of dye seqs, a TSV file of
// radiometries, and write predicted classifications to a TSV file.

#include <fstream>
#include <iomanip>  // for std::setprecision
#include <iostream>
#include <string>

#include "common/dye_seq.h"
#include "common/error_model.h"
#include "common/radiometry.h"
#include "classifiers/fwd_alg_classifier.h"
#include "classifiers/scored_classification.h"

#ifdef _OPENMP
#include <omp.h>
#include "classifiers/omp_fwd_alg_classifier.h"
#else
#include "classifiers/fwd_alg_classifier.h"
#endif

namespace {
using fluoroseq::DistributionType;
using fluoroseq::DyeSeq;
using fluoroseq::ErrorModel;
using fluoroseq::Radiometry;
using fluoroseq::ScoredClassification;
using std::cout;
using std::ifstream;
using std::ofstream;
using std::setprecision;
using std::string;
#ifdef _OPENMP
using fluoroseq::OMPFwdAlgClassifier;
#else
using fluoroseq::FwdAlgClassifier;
#endif
}

int main(int argc, char** argv) {
    if (argc != 4) {
        cout << "wrong number of inputs\n";
        return EXIT_FAILURE;
    }
    char* dye_seqs_filename = argv[1];
    char* radiometries_filename = argv[2];
    char* predictions_filename = argv[3];

    ErrorModel error_model(.06,  // p_edman_failure
                           .05,  // p_detach
                           .05,  // p_bleach
                           .07,  // p_dud
                           DistributionType::LOGNORMAL,
                           1.0,  // mu
                           .16);  // sigma
    cout << "built error model\n";

    ifstream fdye(dye_seqs_filename);
    int num_channels = 0;
    fdye >> num_channels;
    int num_dye_seqs = 0;
    fdye >> num_dye_seqs;
    DyeSeq** dye_seqs = new DyeSeq*[num_dye_seqs];
    for (int i = 0; i < num_dye_seqs; i++) {
        string dye_string;
        fdye >> dye_string;
        int num_peptides;
        fdye >> num_peptides;
        int id;
        fdye >> id;
        dye_seqs[i] = new DyeSeq(num_channels, dye_string, num_peptides, id);
    }
    fdye.close();
    cout << "read from dye seqs file\n";

    ifstream frad(radiometries_filename);
    int num_timesteps = 0;
    frad >> num_timesteps;
    // int num_channels = 0;
    frad >> num_channels;
    int num_radiometries = 0;
    frad >> num_radiometries;
    Radiometry** radiometries = new Radiometry*[num_radiometries];
    for (int i = 0; i < num_radiometries; i++) {
        radiometries[i] = new Radiometry(num_timesteps, num_channels);
        for (int j = 0; j < num_timesteps * num_channels; j++) {
            frad >> radiometries[i]->intensities[j];
        }
    }
    frad.close();
    cout << "read from radiometries file\n";

    #ifdef _OPENMP
    OMPFwdAlgClassifier classifier(omp_get_max_threads(),
                                   num_timesteps,
                                   num_channels,
                                   error_model,
                                   num_dye_seqs,
                                   dye_seqs);
    #else
    FwdAlgClassifier classifier(num_timesteps,
                                num_channels,
                                error_model,
                                num_dye_seqs,
                                dye_seqs);
    #endif
    cout << "constructed classifier\n";
    ScoredClassification** results = classifier.classify(num_radiometries,
                                                         radiometries);
    cout << "got results\n";
    
    ofstream fpred(predictions_filename);
    fpred << "radmat_iz,best_pep_iz,best_pep_score\n";
    for (int i = 0; i < num_radiometries; i++) {
        fpred << i << ",";
        fpred << results[i]->y->id << ",";
        fpred << setprecision(17) << results[i]->adjusted_score() << "\n";
    }
    fpred.close();
    cout << "finished\n";
    return 0;
}
