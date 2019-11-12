Varoom
------

Varoom is a C++ toolkit for processing genomic information.  It contains
components for dealing with common data formats such as FASTA, FASTQ,
SAM, BAM, and VCF.  It also contains components for manipulating HGVS
notation variants, and k-mer tools.

These components are arranged in various ways to produce tools for
addressing bioinformatic analysis problems.  The tools are available in
a Docker image.

KL-BAM
======

KL-BAM is a tool for analysing groups of samples that have been sequenced
with the same protocol and identifying loci where individual samples
contain statistically significant variation from the population. Its
intended use is for identifying low-VAF variants in somatic cancer
samples.

