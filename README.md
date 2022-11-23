# PSMC_boostrap_231122
Generate consensus diploid genome and run PSMC plus bootstrap

Description:
A consensus sequence of the BAM file is generated in fastq format using SAMtools mpileup command with the –C50 option to reduce the effect of reads with excessive mismatches; bcftools view –c to call variants; vcfutils.pl vcf2fq to convert the VCF file of called variants to FASTQ format with further filtering to remove sites with less than a third or more than double the average depth of coverage and Phred quality scores less than 30.

PSMC inference is carried out using the recommended input parameters for human autosomal data, i.e. 25 iterations, with maximum TMRCA (Tmax) = 15, number of atomic time intervals (n) = 64 (following the pattern (1*4 + 25*2 + 1*4 + 1*6), with 100 bootstraps (parameters v9-v11). The plot is scaled to time in years based on a given mutation rate and generation time (parameters v7-v8).

PSMC reference: Li, H., Durbin, R., 2011. Inference of human population history from individual whole-genome sequences. Nature 475, 493-496. https://www.nature.com/articles/nature10231

requires local installation of PSMC tools (e.g., ~/programs/psmc): https://github.com/lh3/psmc, https://anaconda.org/genomedk/psmc

##################
Required programs:
* Samtools
* bcftools
* psmc

Input files:
* bam file (reads assembled to reference genome, preferably repeat-masked)
* reference genome fasta file

Use:
The bash script is written to run with SLURM on the NOAA SEDNA genomic analysis cluster.

The bash script has 11 variable parameters to be modified for each PSMC analysis. Variables v1-8 are specific to each assembly and species. Variables v9-11 are specific to PSMC.

The script creates a sudirectory called "PSMC" within the directory containing the bam file, then generates the diploid consensus genome from the bam file. PSMC analysis is run on the diploid consensus genome, and then bootstrapped 100 times. The PSMC output are plotted from the *out.gp file to generate *out.eps files. Plots can be modified (e.g., change axis ranges or plot colors) by editing the *out.gp files and re-plotting using gnuplot. 
