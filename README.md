# umify
Script for UMI error correction

PCR and sequencing errors are corrected by UMI pattern matching and selecting reads containing valid UMIs. Reads with matching UMIS are corrected to the most abundant read per UMI. To exclude quantitative bias introduced during Multiplex PCR reads are normalized to unique mRNA sequences. Reads shorter than 50 nucleotides are filtered for additional quality control.

>perl umify_normalized_iso.pl sample.fasta

The script creates a folder containing a log file with primer and UMI information as well as a new fasta files with corrected sequences for all reads and primer defined isotypes.
