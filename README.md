# Whatprot

A fluoro-sequencing pipeline by Matthew Beauregard Smith

Professor Edward M. Marcotte's lab.

The University of Texas at Austin

Oden Institute and Institute for Cellular and Molecular Biology

## Prerequisites

You will need a linux (or WSL) system with the following:
* g++
* GNU Make (other make programs will likely work, but we use GNU Make)
* python3

## Building from source

Built using g++.

```bash
$ cd whatprot/cc_code
$ make release
```

This will result in a binary `whatprot/cc_code/bin/release/whatprot`.

## Example: classification from a radiometry file

Here we assume you have a data file with radiometries produced by Erisyon's sigproc code, a .fasta file containing proteins you wish to sequence against, and a .json file containing the parameterization of your data.

### Convert the radiometry file.
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

### Produce a dye-seq file starting with a .fasta file.
To produce a 'dye-seq' file that you will use for classification and/or simulation of training data, run the following int he Python REPL of your choice, from the whatprot/python directory.
```python
from cleave_proteins import cleave_proteins
# - Replace filepaths as appropriate.
# - Replace "name-of-protease" with one of the following: "trypsin", "cyanogen bromide", or
#   "EndoPRO".
# - n is optional. If absent every protein in the .fasta file is cleaved, but if present, n
#   proteins are selected randomly without replacement from the provided file.
cleave_proteins("path/to/fasta/file.fasta",
                "path/to/cleaved/peptides/output/file.csv",
                "name-of-protease"
                n=10)

from dye_seqs_from_peptides import dye_seqs_from_peptides
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

### Generate training data (i.e., dye-tracks) if using kNN or hybrid
If you plan to use the kNN or hybrid classifiers, you will need dye-tracks. These can be
simulated by whatprot.
```bash
# Generate dyetrack samples:
#   -t (or --timesteps) number of timesteps to generate during simulation. Should equal number of Edmans + 1.
#   -g (or --numbenerate) number of reads to simulate for each dye-seq. We recommend setting this to 1000.
#   -P (or --seqparams) path to .json file with the sequencing parameters.
#   -S (or --dyeseqs) path to dye-seq file from previous step to generate dye-tracks based on.
#   -T (or --dyetracks) path to dye-track file to save results to.
$ ./bin/release/whatprot simulate dt -t 10 -g 1000 -P ./path/to/parameters.json -S ./path/to/dye-seqs.tsv -T ./path/to/dye-tracks.tsv
```

### Run classification
To classify data using our recommended parameters with the hybrid classifier, run the following:
```bash
# Classify data using the hybrid classifier
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
$ ./bin/release/whatprot classify hybrid -k 10000 -s 0.5 -H 1000 -p 5 -S ./path/to/dye-seqs.tsv -T ./path/to/dye-tracks.tsv -R ./path/to/radiometries.tsv -Y ./path/to/predictions.tsv
```

For smaller datasets we recommend the HMM classifier. In this case run the following:
```bash
# Classify data using the HMM classifier
#   -p (or --hmmprune) pruning cutoff for HMM (measured in sigma of fluorophore/count
#      combination). This parameter is optional; if omitted, no pruning cutoff will be
#      used.
#   -S (or --dyeseqs) dye-seqs to use as reference for HMM classification.
#   -R (or --radiometries) radiometries to classify.
#   -Y (or --results) output file with a classification id and score for every radiometry.
$ ./bin/release/whatprot classify hmm -p 5 -S ./path/to/dye-seqs.tsv -R ./path/to/radiometries.tsv -Y ./path/to/predictions.tsv
```

The kNN classifier is also available for your use. Although we believe the hybrid and HMM classifiers should be better for all use cases, if you choose to run the kNN classifier you can run it as follows:
```bash
# Classify data using the hybrid classifier
#   -k (or --neighbors) number of neighbors to use for kNN part of hybrid classifier.
#   -s (or --sigma) sigma value for gaussian weighting function for neighbor voting.
#   -T (or --dyetracks) dye-tracks to use as training data for kNN classification.
#   -R (or --radiometries) radiometries to classify.
#   -Y (or --results) output file with a classification id and score for every radiometry.
$ ./bin/release/whatprot classify hybrid -k 10000 -s 0.5 -T ./path/to/dye-tracks.tsv -R ./path/to/radiometries.tsv -Y ./path/to/predictions.tsv
```

### Plot precision recall curves

To plot one PR curve for read-level precision and recall run the following in Python
```python
from pr_curve import pr_curve
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
from pr_curve import pr_curve
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

## Performance prediction with simulated data

If you wish to predict the performance by running on simulated test data, you will do the same as described above, but for your radiometries input you will instead generate the data, as follows:
```bash
# Generate dyetrack samples:
#   -t (or --timesteps) number of timesteps to generate during simulation. Should equal number of Edmans + 1.
#   -g (or --numbenerate) number of reads to simulate total. We recommend setting this to 10000.
#   -P (or --seqparams) path to .json file with the sequencing parameters.
#   -S (or --dyeseqs) path to dye-seq file from previous step to generate dye-tracks based on.
#   -R (or --radiometries) path to radiometries file to save results to.
#   -Y (or --results) path to file to save true-ids of the peptides of the generated radiometries.
$ ./bin/release/whatprot simulate dt -t 10 -g 10000 -P ./path/to/parameters.json -S ./path/to/dye-seqs.tsv -R ./path/to/radiometries.tsv -Y ./path/to/true-ids.tsv
```

## File formats

### dye_seqs.tsv

Tab delimited file.

First two lines are special standalone ints that represent: n_channels and n_dye_seqs

Followed by `n_dye_seqs` lines each of which is tab delimited with three terms:
* a "dye_string"
* the number of peptides that generated this dye_seq.
* the id of this dye_seq.

Note that when the id of a dye_seq is produced from peptides, the dye seq is assigned an id
from one of the peptides. The c++ code doesn't know how to read peptide data
though so may not be super relevant, just helps with analysis to be able
to guess at a particular peptide.

#### dye_string
A "dye_string" is a summary of a labeling.

Example dye_string:
```
..0.1  = A peptide that is labelled in the [2] position with channel 0 and [4] with channel 1
```


Example dye_seqs.tsv with 3 channels, 2 sequences:

```
3
2
..0..1.0	2	1
.220........0	3	2
```
