#!/bin/bash
set -e;
export LC_NUMERIC=C;

# Plot original lattice in PDF
lattice-to-fst --acoustic-scale=1 --lm-scale=1 ark:input.txt ark,t:- | \
    tail -n+2 | fstcompile | \
    fstdraw --portrait --isymbols=symbs.txt --osymbols=symbs.txt | \
    dot -Tpdf > input.pdf

# Create output lattice representing the same language as in input.txt
# by removing the CTC blank output labels.
../lattice-remove-ctc-blank 1 ark:input.txt ark,t:output.txt;

# Plot output lattice in PDF
lattice-to-fst --acoustic-scale=1 --lm-scale=1 ark:output.txt ark,t:- | \
    tail -n+2 | fstcompile | \
    fstdraw --portrait --isymbols=symbs.txt --osymbols=symbs.txt | \
    dot -Tpdf > output.pdf
