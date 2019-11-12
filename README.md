Varoom
======

Varoom is a C++ toolkit for processing genomic information.  It contains
components for dealing with common data formats such as FASTA, FASTQ,
SAM, BAM, and VCF.  It also contains components for manipulating HGVS
notation variants, and k-mer tools.

These components are arranged in various ways to produce tools for
addressing bioinformatic analysis problems.  The tools are available in
a Docker image.

KL-BAM
------

KL-BAM is a tool for analysing groups of samples that have been sequenced
with the same protocol and identifying loci where individual samples
contain statistically significant variation from the population. Its
intended use is for identifying low-VAF variants in somatic cancer
samples.

A typical workflow will look like the following:

```bash
    #!/bin/bash

    # Count the base coverage at each position
    #
    mkdir -p sample-counts
    for fn in $(cat my-sample-filenames.txt)
    do
        klbam count --regions my-genes.bed -o sample-counts/$(basename ${fn} .bam).tsv.gz ${fn}
    done

    # Merge the counts to get a global distribution
    #
    klbam merge -o all-counts.tsv.gz sample-counts/*.tsv.gz

    # Compute the KL-divergence from the global distribution for each sample
    #
    mkdir -p sample-dists
    for fn in sample-counts/*.tsv.gz
    do
        klbam dist --global all-counts.tsv.gz ${fn} sample-dists/$(basename ${fn})
    done

    # Estimate distribution parameters from the divergences
    #
    klbam fit -o all-gamma.tsv.gz sample-dists/*.tsv.gz

    # Compute p-values for each sample/position
    #
    mkdir -p sample-pvals
    for fn in sample-dists/*.tsv.gz
    do
        klbam dist --gamma all-gamma.tsv.gz ${fn} sample-pvals/$(basename ${fn})
    done
```

The resulting TSV files can then be filtered according to the desired p-value cutoff,
minimum coverage, minimum VAF, and so on. We love `awk`, and in `R`, `data.table`,
but use whatever works for you!

