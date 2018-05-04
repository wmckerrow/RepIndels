To use:
1. python setup.py
2. Adjust parameters at top of RepIndels.py
3. python RepIndels.py repeats.fa reads_R1.fastq reads_R2.fastq output.fa n_steps

where
repeats.fa is a multi fasta file with all the reapeat sequences from the family of interest,
reads_R1.fastq and reads_R2.fastq are the two fastq files for the paired end reads that align to the repeat of interest
output.fa is where the predicted TE sequences are written in fasta format
n_steps is the number of EM steps to run
