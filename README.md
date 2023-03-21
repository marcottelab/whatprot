# Whatprot

A fluoro-sequencing pipeline by Matthew Beauregard Smith

Professor Edward M. Marcotte's lab.

The University of Texas at Austin

Oden Institute and Institute for Cellular and Molecular Biology

## Table of Contents

* [Prerequisites](#prerequisites)
* [Building from source](#buildingfromsource)
* [Running classification](#runningclassification)
  * [HMM classification](#hmmclassification)
  * [Hybrid classification](#hybridclassification)
  * [kNN classification](#knnclassification)
  * [Multithreaded performance](#multithreadedperformance)
* [Fitting parameters](#fittingparameters)
* [Filetypes](#filetypes)
  * [Sequencing parameters file](#sequencingparametersfile)
  * [Radiometry files](#radiometryfile)
    * [Radiometry file format](#radiometryfileformat)
    * [Converting from Erisyon's format](#radiometryfromerisyon)
    * [Simulating radiometries](#simulatingradiometries)
  * [Dye-seq files](#dyeseqfiles)
    * [Dye-seq file format](#dyeseqfileformat)
    * [Produce a dye-seq file starting with a .fasta file.](#dyeseqfromfasta)
  * [Dye-track files](#dyetrackfiles)
    * [Dye-track file format](#dyetrackfileformat)
    * [Generate dye-tracks](#generatedyetracks)
* [Plotting results](#plottingresults)
* [Sample runthrough with simulated data](#samplerunthroughsimulated)

## Prerequisites <a name='prerequisites' />

You will need a linux (or WSL) system with the following:
* g++
* GNU Make (other make programs will likely work, but we use GNU Make)
* python3 (for converting files from Erisyon formats to whatprot formats and examining results; not necessary for classification).

## Building from source <a name='buildingfromsource' />

Built using g++.

```bash
$ cd whatprot/cc_code
$ make release
```

This will result in a binary `whatprot/cc_code/bin/release/whatprot`.

## Running classification <a name='runningclassification' />
Classify has three different modes: k-Nearest Neighbors (kNN or just NN), Hidden Markov Models (HMM), or hybrid which combines the two approaches. Each needs a different combination of input files. Below we describe the primary use-cases of these differing classification methods and show examples of how to run them. In the "filetypes" section further down, we describe the various filetypes that you may need as input, as well as how to get them. These often come in two forms; (1) simulated data, useful for validating the efficacy of a labeling strategy for a particular organism, and (2) real data, which generally is taken as output from Erisyon's sigproc pipeline.

### HMM classification <a name='hmmclassification' />
For datasets with smaller numbers of peptides we recommend the HMM classifier. Note that this classifier will have a runtime proportional to the product of your number of peptides and the number of reads you want to analyze; this tends to be unreasonable when your reference database has more than perhaps 1000 peptides, though this may depend on your labeling scheme and other factors. If running the HMM classifier we recommend to run it in the following manner:
```bash
# Classify data using the HMM classifier
#   -P (or --seqparams) path to .json file with parameterization information.
#   -p (or --hmmprune) pruning cutoff for HMM (measured in sigma of fluorophore/count
#      combination). This parameter is optional; if omitted, no pruning cutoff will be
#      used.
#   -S (or --dyeseqs) dye-seqs to use as reference for HMM classification.
#   -R (or --radiometries) radiometries to classify.
#   -Y (or --results) output file with a classification id and score for every radiometry.
$ ./bin/release/whatprot classify hmm -p 5 -P ./path/to/seq-params.json -S ./path/to/dye-seqs.tsv -R ./path/to/radiometries.tsv -Y ./path/to/predictions.csv
```

### hybrid classification <a name='hybridclassification' />
The hybrid classifier will greatly improve runtime performance for larger datasets, with little to no impact on the accuracy of your results. To classify data using our recommended parameters with the hybrid classifier, run something like the following:
```bash
# Classify data using the hybrid classifier
#   -P (or --seqparams) path to .json file with parameterization information.
#   -k (or --neighbors) number of neighbors to use for kNN part of hybrid classifier.
#   -s (or --sigma) sigma value for gaussian weighting function for neighbor voting.
#   -H (or --passthrough) max-cutoff for number of peptides to forward from kNN to HMM
#   -p (or --hmmprune) pruning cutoff for HMM (measured in sigma of fluorophore/count
#      combination). This parameter is optional; if omitted, no pruning cutoff will be
#      used.
#   -S (or --dyeseqs) dye-seqs to use as reference for HMM classification.
#   -T (or --dyetracks) dye-tracks to use as training data for kNN classification.
#   -R (or --radiometries) radiometries to classify.
#   -Y (or --results) output file with a classification id and score for every radiometry.
$ ./bin/release/whatprot classify hybrid -k 10000 -s 0.5 -H 1000 -p 5 -P ./path/to/seq-params.json -S ./path/to/dye-seqs.tsv -T ./path/to/dye-tracks.tsv -R ./path/to/radiometries.tsv -Y ./path/to/predictions.csv
```

### kNN classification <a name='knnclassification' />
The kNN classifier is also available for your use. It will be slightly faster than the hybrid classifier but will give much worse results. Although we believe the hybrid and HMM classifiers should be better for all use cases, if you choose to run the kNN classifier you can run it as follows:
```bash
# Classify data using the kNN classifier
#   -P (or --seqparams) path to .json file with parameterization information.
#   -k (or --neighbors) number of neighbors to use for kNN part of hybrid classifier.
#   -s (or --sigma) sigma value for gaussian weighting function for neighbor voting.
#   -T (or --dyetracks) dye-tracks to use as training data for kNN classification.
#   -R (or --radiometries) radiometries to classify.
#   -Y (or --results) output file with a classification id and score for every radiometry.
$ ./bin/release/whatprot classify nn -k 10000 -s 0.5 -P ./path/to/seq-params.json -T ./path/to/dye-tracks.tsv -R ./path/to/radiometries.tsv -Y ./path/to/predictions.csv
```

### Multithreaded performance <a name='multithreadedperformance' />

Classification is automatically multithreaded with OpenMP. You can change the number of threads by setting the OMP_NUM_THREADS environment variable. You may experience sub-optimal performance as the number of threads increases on many linux systems. This is because the default memory allocators included with many linux systems have very poor performance with multithreaded workloads. This issue can be alleviated by using the LD_PRELOAD environment variable to inject a better performing implementation of malloc into the application. We use jemalloc, but we expect that any malloc implementation designed to deal with memory allocations from large numbers of threads will be roughly equivalent in performance (i.e., tcmalloc, hoard, or ptmalloc2).

## Fitting parameters - using real data to determine the correct parameterization. <a name='fittingparameters' />

You will provide a parameter json file in which you provide starting values for your parameters. Note that:
* mu, sig, and bg_sig will be held constant, and should be fit separately before running fit with whatprot.
* Any probability set to 0 or to 1 will be fixed; this is a limitation inherent in the fit procedure being used.

The simplest fitting procedure might look as follows:
```bash
# Fit data using whatprot.
#   -P (or --seqparams) path to .json file with parameterization information.
#   -L (or --stoppingthreshold) iteration will stop when the difference between
#      iterations is less than this value. Note that the difference between the fit
#      value and the true value may be more than this limit.
#   -x (or --dyeseqstring) a sequence of digits and dots (i.e., period characters)
#      representing the peptide being fit. Dots represent amino acids that can't be
#      labeled, while each digit represents an amino acid which can be labeled by a
#      particular type of fluorophore given your labeling scheme. For example '..0.1'
#      is a peptide that can be labelled in position 2 on channel 0 and position 4 on
#      channel 1.
#   -R (or --radiometries) radiometries to fit.
$ ./bin/release/whatprot fit -P /path/to/seq-params.json -L 0.00001 -x ..0.1 -R /path/to/radiometries.tsv
```

Sometimes a confidence interval is desired. You must then additionally specify both the number of bootstrap rounds you want to run and the size of the confidence interval you wish to estimate. You can do this as follows:
```bash
# Fit data using whatprot.
# See previous example for repeated parameters. Additional parameters are:
#   -b (or --numbootstrap) declares the number of bootstrap rounds to run.
#   -c (or --confidenceinterval) declares the size of the confidence interval. Giving
#      0.9 will give a 90% confidence interval using the percentile method, starting
#      at the 5th percentile and ending at the 95th.
$ ./bin/release/whatprot fit -P /path/to/seq-params.json -L 0.00001 -x ..0.1 -R /path/to/radiometries.tsv -b 200 -c .9
```

If you use bootstrapping to produce a confidence interval, you can also get an output file with a table containing the estimated parameter values from every bootstrap run. You would do this as follows:
```bash
# Fit data using whatprot.
# See previous example for repeated parameters. Additional parameter is:
#   -Y (or --results) for the path to where you want to save the bootstrapping information.
$ ./bin/release/whatprot fit -P /path/to/seq-params.json -L 0.00001 -x ..0.1 -R /path/to/radiometries.tsv -b 200 -c .9 -Y /path/to/results.csv
```

## Filetypes - what they are and how to get them. <a name='filetypes' />

### Sequencing parameters file - contains your parameterization of the sequencing process. <a name='sequencingparametersfile' />

Here is an example file with one channel:

```json
{
  "p_edman_failure": 0.06,
  "p_detach": 0.05,
  "p_initial_block": 0.07,
  "p_cyclic_block": 0.02,
  "channel_models": [
    {
      "p_bleach": 0.05,
      "p_dud": 0.07,
      "bg_sig": 0.00667,
      "mu": 1.0,
      "sig": 0.16
    }
  ]
}
```

It's easy to add more channels, for example to have two do this:

```json
{
  "p_edman_failure": 0.06,
  "p_detach": 0.05,
  "p_initial_block": 0.07,
  "p_cyclic_block": 0.02,
  "channel_models": [
    {
      "p_bleach": 0.05,
      "p_dud": 0.07,
      "bg_sig": 66.7,
      "mu": 10000,
      "sig": 1600
    },
    {
      "p_bleach": 0.05,
      "p_dud": 0.07,
      "bg_sig": 66.7,
      "mu": 10000,
      "sig": 1600
    }
  ]
}
```

### Radiometry file - contains 'reads' you wish to classify. <a name='radiometryfile' />

Will have .tsv (tab-separated values) filetype.

#### Radiometry file format <a name='radiometryfileformat' />

The first three lines are special.
* Fist line is an integer representing the number of timesteps (i.e., number of Edman cycles plus one).
* Second line is an integer representing the number of channels (number of colors or types of fluorophore).
* Third line is an integer representing the number of reads in the file.

The rest of the lines in the file contain read data. Each line is one read. Tab-separated values represent intensity values at each channel and timestep. Note that each line iterates first through channels and then through timesteps, i.e., timestep 0 channel 0 TAB timestep 0 channel 1 TAB timestep 1 channel 0..... This is the opposite arrangement to the organization of data from Erisyon. This is handled by the convert_radiometries python function described below.

Example radiometries.tsv
```
5
2
3
1050.0  120.0  970.0  -30.0 1020.0  10.0  -40.0 20.0  70.0  60.0
10.0  1010.0  30.0  970.0  20.0  1020.0  -40.0 1040.0  50.0  920.0
-30.0 1070.0  130.0  860.0  90.0  1030.0  -60.0 60.0  20.0  -110.0
```

#### Converting from Erisyon's format <a name='radiometryfromerisyon' />

You will need to convert your radiometries from the format used by Erisyon to the format used by whatprot. Run the following int he Python REPL of your choice, from the whatprot/python directory.
```python
from convert_radiometries import convert_radiometries
# - num_channels is the number of colors of fluorophore.
# - num_mocks is the number of mocks in the input file. This many cycles will be removed
#   from the beginning of every read.
# - num_cycles is the number of Edmans plus one (i.e., including the 'pre-edman' cycle
#   which you should always sequence with).
convert_radiometries(num_channels,
                     num_mocks,
                     num_cycles,
                     "path/to/erisyon/radmat/file.npy/or/.tsv",
                     "path/to/whatprot/format/radmat/file.tsv")
```

#### Simulating radiometries <a name='simulatingradiometries' />

If you wish to predict the performance by running on simulated test data, you will do the same as described above, but for your radiometries input you will instead generate the data, as follows:
```bash
# Generate radiometry samples:
#   -t (or --timesteps) number of timesteps to generate during simulation. Should equal number of Edmans + 1.
#   -g (or --numgenerate) number of reads to simulate total. We recommend setting this to 10000. The actual
#      number of reads will be less, because reads of unlabelable peptides or all dud-dyes will be removed
#      from the results -- these would not be seen in real data.
#   -P (or --seqparams) path to .json file with the sequencing parameters.
#   -S (or --dyeseqs) path to dye-seq file from previous step to generate dye-tracks based on.
#   -R (or --radiometries) path to radiometries file to save results to.
#   -Y (or --results) path to file to save true-ids of the peptides of the generated radiometries.
$ ./bin/release/whatprot simulate rad -t 10 -g 10000 -P ./path/to/parameters.json -S ./path/to/dye-seqs.tsv -R ./path/to/radiometries.tsv -Y ./path/to/true-ids.tsv
```

### dye-seq files - contains abstract representation of sequenceable information given a labeling scheme. <a name='dyeseqfiles'/>

Will have .tsv (tab-separated values) filetype.

#### Dye-seq file format <a name='dyeseqfileformat' />

First two lines are special 
* First line is an integer representing the number of channels (number of colors or types of fluorophore).
* Second line is an integer representing the number of dye-seqs in the file.

Followed by `n_dye_seqs` lines each of which is tab delimited with three terms:
* a "dye_string" which is a sequence of digits and dots (i.e., period characters). Dots represent amino acids that can't be labeled, while each digit represents an amino acid which can be labeled by a particular type of fluorophore given your labeling scheme. For example '..0.1' is a peptide that can be labelled in position 2 on channel 0 and position 4 on channel 1.
* the number of peptides that generated this dye_seq.
* the id of this dye_seq, which will be the lowest valued id of the peptides which generate it.

Example dye_seqs.tsv with 3 channels, 2 sequences:

```
3
2
..0..1.0	2	1
.220........0	3	2
```

#### Produce a dye-seq file starting with a .fasta file. <a name='dyeseqfromfasta' />
To produce a 'dye-seq' file that you will use for classification and/or simulation of training data, run the following int he Python REPL of your choice, from the whatprot/python directory.
```python
from cleave_proteins import cleave_proteins
from dye_seqs_from_peptides import dye_seqs_from_peptides

# - Replace filepaths as appropriate.
# - Replace "name-of-protease" with one of the following: "trypsin", "cyanogen bromide", or
#   "EndoPRO".
# - n is optional. If absent every protein in the .fasta file is cleaved, but if present, n
#   proteins are selected randomly without replacement from the provided file.
cleave_proteins("path/to/fasta/file.fasta",
                "path/to/cleaved/peptides/output/file.csv",
                "name-of-protease"
                n=10)

# - Replace filepaths as appropriate.
# - label_set should be a list of strings, each string containing a set of characters
#   (typically but not always just one character). Each index of the list represents a
#   channel or color of fluorophore, and each letter represents an amino acid code.
#   For example, if we use the label_set ['DE','C','Y'], then we will label aspartic
#   acid (D) and glutamic acid (E) on channel 0, cysteine (C) on channel 1, and tyrosine
#   (Y) on channel 2.
dye_seqs_from_peptides("path/to/peptides/file/from/previous/step.csv",
                       label_set,
                       "path/to/dye-seqs/output/file.tsv")
```

### Dye-track files - contains training data for kNN or hybrid classifiers <a name='dyetrackfiles' />

You will always simulate these, and you will need a dye-seq to do so.

#### Dye-track file format <a name='dyetrackfileformat' />

First three lines are special
* First line is number of timesteps
* Second line is number of channels
* Third line is number of dye-tracks.

After that each line represents one dye-track. Values are tab delimited.
* First num timesteps times num channels values are integer values representing the number of remaining fluorophores of that channel at that timestep. This is organized in the same manner as the radiometries file; first by channel, then by timestep.
* After that there is a number representing the number of distinct dye-seqs which produced this dye-track during simulation.
* What follows is a list of that length times three; each following triple of values (still tab separated), is first the ID of the dye-seq, then the number of peptides mapping to that dye seq, and thirdly the number of 'hits' or simulated reads from that dye-seq produced the dye-track in this row.

Example dye-track file:
```
5
2
3
1 0 1 0 1 0 1 0 0 0 1 171 1 3
1 2 0 2 0 2 0 1 0 0 3 171 1 5 235 3 8 113 1 1
0 2 0 2 0 1 0 1 0 1 2 116 1 1 115 1 1
```

#### Generate dye-tracks <a name='generatedyetracks' />
```bash
# Generate dyetrack samples:
#   -t (or --timesteps) number of timesteps to generate during simulation. Should equal number of Edmans + 1.
#   -g (or --numgenerate) number of reads to simulate for each peptide (more generated for dye-seqs mapping to
#      multiple peptides). We recommend setting this to 1000.
#   -P (or --seqparams) path to .json file with the sequencing parameters.
#   -S (or --dyeseqs) path to dye-seq file from previous step to generate dye-tracks based on.
#   -T (or --dyetracks) path to dye-track file to save results to.
$ ./bin/release/whatprot simulate dt -t 10 -g 1000 -P ./path/to/parameters.json -S ./path/to/dye-seqs.tsv -T ./path/to/dye-tracks.tsv
```

## Plotting results <a name='plottingresults' />

To plot one PR curve for read-level precision and recall run the following in Python
```python
from pr_curve import plot_pr_curve
import matplotlib.pyplot as plt

# - The predictions file is your classification results.
# - The true values is a file with one line for each entry in your predictions file
#   set to the value you expect as a result. For verifying the algorithm on simulated
#   data (see below), this is generated when you generate your test data.
# - The dye-seqs file is the set of dye-seqs you are computing on. This is necessary
#   because the number of peptides for each dye-seq is needed to properly weight the
#   results.
# - their is an optional "directory" parameter for convenience. If given it will be
#   prepended to all of your filepaths.
plot_pr_curve("path/to/predictions.tsv",
              "path/to/true/values.tsv",
              "path/to/dye-seqs/file.tsv")
```

To plot one PR curve for each of read-level, peptide-level, and protein-level precision
and recall run the following in Python.
```python
from pr_curve import plot_pr_curve
import matplotlib.pyplot as plt

# - First three variables are as before.
# - full_peps_file and lim_peps_file should always either both be specified or neither.
# - full_peps_file is the file with cleaved peptides you used to generate dye-seqs. You 
#   MUST use the same file or this may not work correctly.
# - lim_peps_file is the file with only the peptides you expect to see. This needs to be
#   based on the same set of proteins as your full_peps file (i.e., protein ID for the
#   same peptides should match in the two files. This is the second column in the .csv).
plot_pr_curve("path/to/predictions.tsv",
              "path/to/true/values.tsv",
              "path/to/dye-seqs/file.tsv",
              full_peps_file="path/to/full/peptides/file/including/decoys.csv"
              lim_peps_file="path/to/target/or/true-set/peptides.csv")
```

To plot multiple PR curves together, use instead the 'plot_pr_curves()' function in pr_curves.py (note the extra 's' at the end of the function name). The parameter ordering is the same. You must then provide a list of prediction files instead of just one. You may optionally provide a list of true-values files, a list of dye-seqs files, or even lists of full peptides files or limited (true-set) peptide files. For each of these variables, if one value is given it is used for every predictions file specified, and if you instead provide a list then the values are collated.

## Sample runthrough with simulated data <a name='samplerunthroughsimulated' />

First we make a directory to put our intermediate files and results into, then go into the python subdirectory and start python.
```bash
$ mkdir temp
$ cd python
$ python
```

In python, we will cleave a random set of 100 proteins, and then apply a labeling scheme, in order to get a dye-seqs file.
```python
from cleave_proteins import cleave_proteins
from dye_seqs_from_peptides import dye_seqs_from_peptides
cleave_proteins("../examples/UP000005640_9606.fasta",
                "../temp/peptides.tsv",
                "trypsin",
                n=100)
dye_seqs_from_peptides("../temp/peptides.tsv",
                       ['DE','C','Y'],
                       "../temp/dye-seqs.tsv")
```

After exiting python, we go into the cc_code directory and simulate radiometries (simulate rad), simulate dye-tracks (simulate dt), classify using the hybrid model (classify hybrid) and then go back to the python directory and start python (this section assumes you have already built the application - if not see [Building from source](#buildingfromsource)).
```bash
$ cd ../cc_code
$ ./bin/release/whatprot simulate rad -t 10 -g 10000 -P ../examples/seqparams_atto647n_x3.json -S ../temp/dye-seqs.tsv -R ../temp/radiometries.tsv -Y ../temp/true-ids.tsv
$ ./bin/release/whatprot simulate dt -t 10 -g 1000 -P ../examples/seqparams_atto647n_x3.json -S ../temp/dye-seqs.tsv -T ../temp/dye-tracks.tsv
$ ./bin/release/whatprot classify hybrid -k 10000 -s 0.5 -H 1000 -p 5 -P ../examples/seqparams_atto647n_x3.json -S ../temp/dye-seqs.tsv -T ../temp/dye-tracks.tsv -R ../temp/radiometries.tsv -Y ../temp/predictions.csv
$ cd ../python
$ python
```

Back in python, we plot our results.
```python
from pr_curve import plot_pr_curve
import matplotlib.pyplot as plt
plot_pr_curve("../temp/predictions.csv",
              "../temp/true-ids.tsv",
              "../temp/dye-seqs.tsv")
plt.show()
```

What you get should look something like this:
![example](https://user-images.githubusercontent.com/3892206/193147900-16c8ec47-5837-4bdc-bb7a-69b59cac09c4.png)

Note that you may not get exactly the same result, as there is randomness involved. Nevertheless your result should be similar.
