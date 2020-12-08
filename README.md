# detect_recombination
This repository has the scripts used by gabaldonlab to detect recombination events in hybrid genomes.

-detect_recombination_PHASEDref.py can only be used when the paired-end reads of a hybrid are simultaneously mapped to both parental haplotypes.

-detect_recombination_PHASING.py can only be used when the paired-end reads of a hybrid are mapped on a single haplotype, which corresponds to a parental lineage. This script requires HapCUT2 to phase the heterozygous variants, and relies on these results to identify recombinant regions (change in genotype from 0|1 to 1|0).
