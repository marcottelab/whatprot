# Whatprot

A fluoro-sequencing pipeline by Matthew Beauregard Smith

Professor Edward M. Marcotte's lab.

The University of Texas at Austin

Oden Institute and Institute for Cellular and Molecular Biology


# Building from source

Built using g++.

```bash
$ cd whatprot/cc_code
$ make
```

This will result in a binary `whatprot/cc_code/bin/debug/whatprot`.

# Examples

Assuming you have a PATH to `whatprot/cc_code/bin/debug`

```bash
# Generate dyetrack samples:
#   10 Edman cycles
#   100 samples per peptide
#   using the TSV file of dye_seqs in ./example_files/dye_seqs.tsv
#   writing output to ./dye_tracks.tsv
$ whatprot simulate dt 10 100 ./example_files/dye_seqs.tsv ./dye_tracks.tsv

```

# File formats

## dye_seqs.tsv

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

### dye_string
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