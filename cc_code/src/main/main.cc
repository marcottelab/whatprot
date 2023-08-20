/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Standard C++ library headers:
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

// External headers:
#include "cxxopts.hpp"

// Local project headers:
#include "main/cmd-line-out.h"
#include "main/run-classify-hmm.h"
#include "main/run-classify-hybrid.h"
#include "main/run-classify-nn.h"
#include "main/run-fit.h"
#include "main/run-score-hmm.h"
#include "main/run-simulate-dt.h"
#include "main/run-simulate-rad.h"

namespace {
using cxxopts::Options;
using cxxopts::ParseResult;
using cxxopts::PositionalList;
using cxxopts::value;
using std::allocator;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using whatprot::print_bad_inputs;
using whatprot::print_invalid_command;
using whatprot::print_omp_info;
using whatprot::run_classify_hmm;
using whatprot::run_classify_hybrid;
using whatprot::run_classify_nn;
using whatprot::run_fit;
using whatprot::run_score_hmm;
using whatprot::run_simulate_dt;
using whatprot::run_simulate_rad;
}  // namespace

int main(int argc, char** argv) {
    Options options("whatprot",
                    "whatprot is a program for analyzing protein "
                    "fluorosequencing data.");

    // clang-format off
    options.add_options()
        ("h,help", "Print usage\n")
        ("b,numbootstrap",
            "Only for parameter fitting, and optional. If specified, indicates "
            "number of bootstrapping rounds to perform to get confidence "
            "intervals for parameters. If specified you must also specify "
            "--confidenceinterval (shorthand -c).\n",
            value<int>())
        ("c,confidenceinterval",
            "Only for parameter fitting, and optional. If specified, indicates "
            "the desired confidence interval size. If specified you must also "
            "specify --numbootstrap (shorthand -b).\n",
            value<double>())
        ("g,numgenerate",
            "Only for simulation, and required. Number of dye-tracks or "
            "radiometries to generate. For simulate rad, this is the actual "
            "number. For simulate dt, this is the number to generate PER "
            "PEPTIDE. Note also that the real number written to the output "
            "file may be less than this. When results would not be visible in "
            "a real run due to failure to attach functioning fluorophores "
            "prior to sequencing, the result is omitted in the output file.\n",
            value<int>())
        ("k,neighbors",
            "Only for nn or hybrid classification, and required. Number of "
            "neighbors to use for kNN classification\n",
            value<int>())
        ("p,hmmprune",
            "Only for hmm or hybrid classification, and NOT required. Defines "
            "a multiplier on sigma to use when pruning an HMM for greater "
            "efficiency. Higher values imply less pruning.\n",
            value<double>())
        ("s,sigma",
            "Only for nn or hybrid classification, and required. Sigma to use "
            "for the Gaussian kernel used to weight votes in kNN "
            "classification.\n",
            value<double>())
        ("t,timesteps",
            "Only for simulation, and required. Number of timesteps to "
            "generate during simulation.\n",
            value<int>())
        ("x,dyeseqstring",
            "Only for fit, and required. Sequence of numerals and dots (.) "
            "indicating the true sequence of fluorophores on the peptide in "
            "the ideal case. Each character represents one amino acid, as if "
            "the right-most character is anchored to the flow-cell. The "
            "numerals indicate which channel, while a '.' indicates an amino "
            "acid with no fluorophore. For example \"..0.1\" is a peptide with "
            "a fluorophore on channel 0 at position 3 and on channel 1 at "
            "position 5.\n",
            value<string>())
        ("H,passthrough",
            "Only for hybrid classification, and required. Number of peptide "
            "candidates to pass through from kNN to HMM.\n",
            value<int>())
        ("L,stoppingthreshold",
            "Only for fit, and required. Threshold of change in the norm to "
            "stop the fitting iteration.\n",
            value<double>())
        ("M,maxruntime",
            "Only for fit, but NOT required. Sets a time limit in minutes, "
            "after which runtime stops, even if the --stoppingthreshold (-L) "
            "has not been reached. Note that this runtime cut-off is "
            "approximate and some parts of the code are exempt, and users "
            "should provide a large margin of error if runtime constraints are "
            "tight. Additionally, if used in combination with bootstrapping "
            "(when --confidenceinterval and --numbootstrap are set), this "
            "runtime limit is used for each bootstrap round and again for the "
            "original data to get a best fit. In that case, the true runtime "
            "may run as high as M * (1 + (1 + b) / t), where M is the value of "
            "--maxruntime, b is the value of --numbootstrap and t is the "
            "number of threads, which is dependent on your specific computing "
            "set-up.\n",
            value<int>())
        ("P,seqparams",
            "Needed by all code paths except for kNN classification. Provides"
            "necessary modeling params for HMM classification (and by extension"
            "hybrid classification), simulation parameters, and starting"
            "estimates for fitting procedures.\n",
            value<string>())
        ("R,radiometries",
            "Only for classification, fit, or rad simulation, and required. "
            "Name of file to read radiometries from (classification or fit) or "
            "name of file to write them to (rad simulation).\n",
            value<string>())
        ("S,dyeseqs",
            "Only for hmm or hybrid classification, and required. Name of file "
            "to read dye-seqs from.\n",
            value<string>())
        ("T,dyetracks",
            "Only for nn or hybrid classification, or dt simulation, and "
            "required. Name of file to read dye-tracks from (classification) "
            "or name of file to write them to (dt simulation).\n",
            value<string>())
        ("Y,results",
            "Required for classify and simulate rad. Then it is the name of "
            "the file to write correct or predicted peptide ids (y values) to. "
            "Optional for fit, and then only permitted if --numbootstrap "
            "(shorthand -b) and --confidenceinterval (shorthand -c) are "
            "specified. Then it is the name of the file to write verbose "
            "information about parameter fit results for each bootstrap run.\n",
            value<string>());
    // clang-format on
    options.custom_help(
            "[MODE] [VARIANT] [OPTS...]\n\n"
            "  MODE is one of classify, fit, or simulate. Other parameter\n"
            "  requirements depend on these.\n"
            "  \n"
            "  For MODE classify, you must define a VARIANT as one of hmm,\n"
            "  hybrid, or nn. Your data will then be classified using the\n"
            "  model defined by VARIANT. These VARIANT possibilities require\n"
            "  specific parameters.\n"
            "  \n"
            "    For VARIANT hmm, you must define --seqparams, --dyeseqs,\n"
            "    --radiometries, and --results. Option --hmmprune is also\n"
            "    permitted.\n"
            "    \n"
            "    For VARIANT hybrid, you must define --seqparams,\n"
            "    --neighbors, --sigma, --passthrough, --dyeseqs, --dyetracks,\n"
            "    --radiometries, and --results. Option --hmmprune is also\n"
            "    permitted.\n"
            "    \n"
            "    For VARIANT nn, you must define --seqparams, --neighbors,\n"
            "    --sigma, --dyetracks, --radiometries, and --results.\n"
            "    \n"
            "  For MODE fit, you must NOT define a VARIANT, and you MUST\n"
            "  define --seqparams, --stoppingthreshold, --dyeseqstring, and\n"
            "  --radiometries. You may specify BOTH --numbootstrap and\n"
            "  --confidenceinterval in order to use bootstrapping to produce\n"
            "  a confidence interval. If you do so (and only if), you may\n"
            "  also specify --results to get complete results for your run.\n"
            "  \n"
            "  For MODE simulate, you must define a VARIANT as either dt or\n"
            "  rad. Your data will then be simulated as dye-tracks or\n"
            "  radiometries depending on your choice. These VARIANT\n"
            "  possibilities require specific parameters.\n"
            "  \n"
            "    For VARIANT dt, you must define --seqparams, --timesteps,\n"
            "    --numgenerate, --dyeseqs, and --dyetracks.\n"
            "    \n"
            "    For VARIANT rad, you must define --seqparams, --timesteps,\n"
            "    --numgenerate, --dyeseqs, --radiometries, and --results\n"
            "    \n");

    // Parse options.
    options.allow_unrecognised_options();
    ParseResult parsed_opts = options.parse(argc, argv);
    if (parsed_opts.count("help")) {
        cout << options.help() << endl;
        return 0;
    }
    const vector<string, allocator<string>>& positional_args =
            parsed_opts.unmatched();

    // Retrieve options, keeping track of which ones were or weren't set.
    unsigned int num_optional_args = 0;
    bool has_b = false;
    int b = -1;
    if (parsed_opts.count("numbootstrap")) {
        has_b = true;
        num_optional_args++;
        b = parsed_opts["numbootstrap"].as<int>();
    }
    bool has_c = false;
    double c = 0.0;
    if (parsed_opts.count("confidenceinterval")) {
        has_c = true;
        num_optional_args++;
        c = parsed_opts["confidenceinterval"].as<double>();
    }
    bool has_g = false;
    int g = -1;
    if (parsed_opts.count("numgenerate")) {
        has_g = true;
        num_optional_args++;
        g = parsed_opts["numgenerate"].as<int>();
    }
    bool has_k = false;
    int k = -1;
    if (parsed_opts.count("neighbors")) {
        has_k = true;
        num_optional_args++;
        k = parsed_opts["neighbors"].as<int>();
    }
    bool has_p = false;
    double p = std::numeric_limits<double>::max();
    if (parsed_opts.count("hmmprune")) {
        has_p = true;
        num_optional_args++;
        p = parsed_opts["hmmprune"].as<double>();
    }
    bool has_s = false;
    double s = 0.0;
    if (parsed_opts.count("sigma")) {
        has_s = true;
        num_optional_args++;
        s = parsed_opts["sigma"].as<double>();
    }
    bool has_t = false;
    int t = -1;
    if (parsed_opts.count("timesteps")) {
        has_t = true;
        num_optional_args++;
        t = parsed_opts["timesteps"].as<int>();
    }
    bool has_x = false;
    string x("");
    if (parsed_opts.count("dyeseqstring")) {
        has_x = true;
        num_optional_args++;
        x = parsed_opts["dyeseqstring"].as<string>();
    }
    bool has_H = false;
    int H = -1;
    if (parsed_opts.count("passthrough")) {
        has_H = true;
        num_optional_args++;
        H = parsed_opts["passthrough"].as<int>();
    }
    bool has_L = false;
    double L = 0.0;
    if (parsed_opts.count("stoppingthreshold")) {
        has_L = true;
        num_optional_args++;
        L = parsed_opts["stoppingthreshold"].as<double>();
    }
    bool has_M = false;
    int M = std::numeric_limits<int>::max();
    if (parsed_opts.count("maxruntime")) {
        has_M = true;
        num_optional_args++;
        M = parsed_opts["maxruntime"].as<int>();
    }
    bool has_P = false;
    string P("");
    if (parsed_opts.count("seqparams")) {
        has_P = true;
        num_optional_args++;
        P = parsed_opts["seqparams"].as<string>();
    }
    bool has_R = false;
    string R("");
    if (parsed_opts.count("radiometries")) {
        has_R = true;
        num_optional_args++;
        R = parsed_opts["radiometries"].as<string>();
    }
    bool has_S = false;
    string S("");
    if (parsed_opts.count("dyeseqs")) {
        has_S = true;
        num_optional_args++;
        S = parsed_opts["dyeseqs"].as<string>();
    }
    bool has_T = false;
    string T("");
    if (parsed_opts.count("dyetracks")) {
        has_T = true;
        num_optional_args++;
        T = parsed_opts["dyetracks"].as<string>();
    }
    bool has_Y = false;
    string Y("");
    if (parsed_opts.count("results")) {
        has_Y = true;
        num_optional_args++;
        Y = parsed_opts["results"].as<string>();
    }

    // Now we determine which variant of the program is being run, and act
    // appropriately based on the options.
    if (positional_args.size() < 1) {
        cout << endl << "INCORRECT USAGE" << endl << endl;
        cout << options.help() << endl;
        return 1;
    }
    if (0 == positional_args[0].compare("classify")) {
        if (positional_args.size() != 2) {
            cout << endl << "INCORRECT USAGE" << endl << endl;
            cout << options.help() << endl;
            return 1;
        }
        if (0 == positional_args[1].compare("hmm")) {
            // Special handling for p since it is optional for classify hmm.
            if (has_p) {
                num_optional_args--;
            }
            if (num_optional_args != 4 || !has_P || !has_S || !has_R
                || !has_Y) {
                cout << endl << "INCORRECT USAGE" << endl << endl;
                cout << options.help() << endl;
                return 1;
            }
            print_omp_info();
            run_classify_hmm(P, p, S, R, Y);
            return 0;
        }
        if (0 == positional_args[1].compare("hybrid")) {
            // Special handling for p since it is optional for classify hybrid.
            if (has_p) {
                num_optional_args--;
            }
            if (num_optional_args != 8 || !has_P || !has_k || !has_s || !has_H
                || !has_S || !has_T || !has_R || !has_Y) {
                cout << endl << "INCORRECT USAGE" << endl << endl;
                cout << options.help() << endl;
                return 1;
            }
            print_omp_info();
            run_classify_hybrid(P, k, s, H, p, S, T, R, Y);
            return 0;
        }
        if (0 == positional_args[1].compare("nn")) {
            if (num_optional_args != 6 || !has_P || !has_k || !has_s || !has_T
                || !has_R || !has_Y) {
                cout << endl << "INCORRECT USAGE" << endl << endl;
                cout << options.help() << endl;
                return 1;
            }
            print_omp_info();
            run_classify_nn(P, k, s, T, R, Y);
            return 0;
        }
        cout << endl << "INCORRECT USAGE" << endl << endl;
        cout << options.help() << endl;
        return 1;
    }
    if (0 == positional_args[0].compare("fit")) {
        // Special handling for b and c. Must have both or neither.
        if (has_b ^ has_c) {
            cout << endl << "INCORRECT USAGE" << endl << endl;
            cout << options.help() << endl;
            return 1;
        }
        // Special handling for b and c. These args are optional.
        if (has_b & has_c) {
            num_optional_args -= 2;
            // If we have b and c, then Y is permitted, and optional.
            if (has_Y) {
                num_optional_args--;
            }
        }
        // Special handling for M since it is optional.
        if (has_M) {
            num_optional_args--;
        }
        if (positional_args.size() != 1 || num_optional_args != 4 || !has_P
            || !has_L || !has_x || !has_R) {
            cout << endl << "INCORRECT USAGE" << endl << endl;
            cout << options.help() << endl;
            return 1;
        }
        print_omp_info();
        // Convert M from more human-readable minutes as integer, to more
        // machine readable seconds as double.
        run_fit(L, 60.0 * (double)M, x, P, R, b, c, Y);
        return 0;
    }
    if (0 == positional_args[0].compare("score")) {
        if (0 == positional_args[1].compare("hmm")) {
            // Special handling for p since it is optional for classify hmm.
            if (has_p) {
                num_optional_args--;
            }
            if (num_optional_args != 3 || !has_S || !has_R || !has_Y) {
                cout << endl << "INCORRECT USAGE" << endl << endl;
                cout << options.help() << endl;
                return 1;
            }
            print_omp_info();
            run_score_hmm(p, S, R, Y);
            return 0;
        }
        cout << endl << "INCORRECT USAGE" << endl << endl;
        cout << options.help() << endl;
        return 1;
    }
    if (0 == positional_args[0].compare("simulate")) {
        if (positional_args.size() != 2) {
            cout << endl << "INCORRECT USAGE" << endl << endl;
            cout << options.help() << endl;
            return 1;
        }
        if (0 == positional_args[1].compare("dt")) {
            if (num_optional_args != 5 || !has_P || !has_t || !has_g || !has_S
                || !has_T) {
                cout << endl << "INCORRECT USAGE" << endl << endl;
                cout << options.help() << endl;
                return 1;
            }
            run_simulate_dt(t, g, P, S, T);
            return 0;
        }
        if (0 == positional_args[1].compare("rad")) {
            if (num_optional_args != 6 || !has_P || !has_t || !has_g || !has_S
                || !has_R || !has_Y) {
                cout << endl << "INCORRECT USAGE" << endl << endl;
                cout << options.help() << endl;
                return 1;
            }
            run_simulate_rad(t, g, P, S, R, Y);
            return 0;
        }
        cout << endl << "INCORRECT USAGE" << endl << endl;
        cout << options.help() << endl;
        return 1;
    }
    cout << endl << "INCORRECT USAGE" << endl << endl;
    cout << options.help() << endl;
    return 1;
}
