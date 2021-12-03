# umify
Script for UMI error correction

Dependencies: Perl5

PCR and sequencing errors are corrected by UMI pattern matching and selecting reads containing valid UMIs. Reads with matching UMIs are corrected to the most abundant read per UMI. To exclude quantitative bias introduced during Multiplex PCR reads are normalized to unique mRNA sequences. Reads shorter than 50 nucleotides are filtered for additional quality control.

## Installation



## command line usage

Please note the script only works with primers described here:
"method paper star protocol"


>perl umify_normalized_iso.pl sample.fasta


Output: The script creates a folder containing a log file with primer and UMI information as well as new fasta files with corrected sequences for all reads and primer defined isotypes.






## Citing umify

If you are using umify in any published work, please cite:

"method paper star protocol"
