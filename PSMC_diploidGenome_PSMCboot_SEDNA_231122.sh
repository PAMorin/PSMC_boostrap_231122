#! /bin/bash

#SBATCH --job-name=PSMC
#SBATCH -e PSMC_%j.e.txt
#SBATCH -o PSMC_%j.log 
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  # (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 20    # <number of cores to ask for>
#SBATCH --mem=80G  
#SBATCH -t 3-0   # walltime in mins or mins:secs, hrs:mins:secs or days-hours

######################
# Phil Morin, 23 Nov. 2022
# A consensus sequence of the BAM file is generated in fastq format using SAMtools mpileup command with the –C50 option to reduce the effect of reads with excessive mismatches; bcftools view –c to call variants; vcfutils.pl vcf2fq to convert the VCF file of called variants to FASTQ format with further filtering to remove sites with less than a third or more than double the average depth of coverage and Phred quality scores less than 30
# PSMC inference is carried out using the recommended input parameters for human autosomal data, i.e. 25 iterations, with maximum TMRCA (Tmax) = 15, number of atomic time intervals (n) = 64 (following the pattern (1*4 + 25*2 + 1*4 + 1*6), with 100 bootstraps. The plot is scaled to time in years based on a given mutation rate and generation time.
# PSMC reference: Li, H., Durbin, R., 2011. Inference of human population history from individual whole-genome sequences. Nature 475, 493-496. https://www.nature.com/articles/nature10231
# requires local installation of PSMC tools (e.g., ~/programs/psmc): https://github.com/lh3/psmc, https://anaconda.org/genomedk/psmc

#programs
module load bio/samtools/1.11
module load bio/bcftools/1.11
module load tools/gnuplot/5.4.1
PSMC=~/programs/psmc

set -eux

##############################################################
# variable parameters (samples, directories, files)
v1=z0076728										#sample ID
v2=/home/pmorin/projects/Bbai/Bbai_map2Mden_Sept22	#BAM directory
v3=z0076728_merged_dedup_noRepeats.bam				#BAM file
v4=30												#average depth of coverage
v5=/home/pmorin/Ref_genomes/Mden/mMesDen1_NCBI_Ref_genome # reference directory
v6=GCA_025265405.1_mMesDen1_primary_haplotype_genomic.fna #reference file
v7=3.47e-9								#mutation rate (substitutions/site/generation)
v8=20									#generation time (e.g., from Taylor et al. 2007, Table 1, T(r=0)=generation length under pre-disturbance conditions. )
v9=15									# default from Li and Durbin is 15
v10="4+25*2+4+6"						# default from Li and Durbin is 4+25*2+4+6
v11=/scratch/pmorin/temp				#for bootstrap files (temp)
##############################################################

cd ${v2}

SAMPLE=${v1}
THREADS=50

#files
REFDIR=${v5}
REFERENCE=${REFDIR}/${v6}

BAMDIR=${v2}
BAM=${BAMDIR}/${v3} # repeat masked assembly

mkdir -p ${BAMDIR}/PSMC
OUTDIR=${BAMDIR}/PSMC # make PSMC folder in bamdir first
DIPLOID=${SAMPLE}_diploid.fq.gz

# diploid genome parameters
MEANDEPTH=${v4} # has to be an integer for mindepth/maxdepth calc below.
let MIN=MEANDEPTH/3  # 1/3x average coverage
let MAX=MEANDEPTH*2  # 2x average coverage

# PSMC parameters
TEMP_DIR=${v11}
MUTRATE=${v7}  #vaquita estimated rate = 5.83E-9/site/gen (11.9yr/gen*5.83E-9=6.94E-8 µ/site/yr). Mutation rate for Ziphiids with ~20yr/gen = 3.47E-9 µ/site/gen.
MUT="`echo $v7 | sed -e 's/-//g'`"   #3.47e9	# for printing mutation rate on plot and in file name
GEN=${v8}
t=${v9}
PSMC_INT=${v10}
#########################################################################
#########################################################################
#########################################################################

### Create a diploid consensus fastq file from the bam file.
samtools mpileup -uf ${REFERENCE} ${BAM} | bcftools call -c --threads ${THREADS} | vcfutils.pl vcf2fq -d${MIN} -D${MAX} > ${OUTDIR}/${DIPLOID}

#########################################################################
# PSMC
ID=${SAMPLE}

# generate psmcfa file from diploid genome file
${PSMC}/utils/fq2psmcfa -q20 ${OUTDIR}/${DIPLOID} > ${OUTDIR}/${ID}_diploid.psmcfa

#generate the psmc file using the default settings for humans (-N25, -t20, -r5).
${PSMC}/psmc -N25 -t${t} -r5 -p ${PSMC_INT} -o ${OUTDIR}/${ID}_${MUT}_t${t}.psmc ${OUTDIR}/${ID}_diploid.psmcfa

#Make psmc plots and adapt the scaling using this psmc file.
nice ${PSMC}/utils/psmc_plot.pl -u ${MUTRATE} -g ${GEN} -RM ${ID}"_"${MUT} ${OUTDIR}/${ID}_${MUT}_t${t}_psmc.out ${OUTDIR}/${ID}_${MUT}_t${t}.psmc

######################
# bootstrap PSMC

# split the PSMC file
${PSMC}/utils/splitfa ${OUTDIR}/${ID}_diploid.psmcfa > ${OUTDIR}/split_${ID}_diploid.psmcfa

# PSMC bootstrap, multithread = 12 (-P 12)
seq 100 | xargs -P ${THREADS} -i ${PSMC}/psmc -N25 -t${t} -r5 -b -p ${PSMC_INT} -o ${TEMP_DIR}/${ID}_round-{}.psmc ${OUTDIR}/split_${ID}_diploid.psmcfa | sh

# cat original.psmc round-*.psmc > combined.psmc
# within each folder containing original sample psmc and bootstrap files
cat ${OUTDIR}/${ID}_${MUT}_t${t}.psmc ${TEMP_DIR}/${ID}_round-*.psmc > ${OUTDIR}/merged_${ID}_t${t}_boot.psmc

# And then plot it: (-RM -> Id placed on plot).
${PSMC}/utils/psmc_plot.pl -u ${MUTRATE} -g ${GEN} -RM "" ${OUTDIR}/merged_${ID}_t${t}_boot.out ${OUTDIR}/merged_${ID}_t${t}_boot.psmc

cp ${BAMDIR}/*PSMCboot_SEDNA.sh ${OUTDIR}

########################################################################################
# for PSMC parameter info, see https://github.com/lh3/psmc
# the `-p' option specifies that there are 64 atomic time intervals and 28 (=1+25+1+1)
# free interval parameters. The first parameter spans the first 4 atomic time
# intervals, each of the next 25 parameters spans 2 intervals, the 27th spans 4
# intervals and the last parameter spans the last 6 time intervals. The `-p' and
# `-t' options are manually chosen such that after 20 rounds of iterations, at
# least ~10 recombinations are inferred to occur in the intervals each parameter
# spans.

### use the custom script provide in the utils folder with PSMC to convert the fastq file to a pseudo-fasta file. It basically codes regions in 100bp bins as being heterozygote or homozygote.

