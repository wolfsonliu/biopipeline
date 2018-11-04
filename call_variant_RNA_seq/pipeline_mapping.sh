#! /bin/bash
#SBATCH --job-name=P_star
#SBATCH --output=$HOME/log/%x.%j_%A_%a.%N.array.out
#SBATCH --error=$HOME/log/%x.%j_%A_%a.%N.array.err
#SBATCH --partition=C032M0512G
#SBATCH --qos=low
#SBATCH --get-user-env
#SBATCH -n 4
#SBATCH --cpu-freq=high
#SBATCH -A hpc
#SBATCH --mail-type=end
#SBATCH --mail-user=
#SBATCH --time=120:00:00

thread=4

basedir=${HOME}/Project
codedir=${basedir}/code
datadir=${basedir}/data
fqdir=${datadir}/clean
mapdir=${datadir}/bam
starrefdir=${datadir}/star
starfa=${starrefdir}/genome.fa

thefa=${homedir}/Reference/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa

Picard="java -jar ${HOME}/.local/bin/picard.jar"
####################

         # 0      1       2      3      4           5             6        7        8         9     10      11     12       13
labels=(FANCC_1 FANCC_2 Mock_1 Mock_2 Untreated_1 Untreated_2 PPIB111_1 PPIB111_2 PPIB111_3 PPIB_1 PPIB_2 PPIB_3 Random_1 Random_2)

label=${labels[${SLURM_ARRAY_TASK_ID}]}

fq1=${label}.1.fq.gz
fq2=${label}.2.fq.gz

if [ ! -e ${mapdir}/notmap ]; then
    mkdir -p ${mapdir}/notmap
fi

if [ ! -e ${mapdir}/log ]; then
    mkdir -p ${mapdir}/log
fi

####################
# star

if [ ! -e ${mapdir}/${label}/star1 ]; then
    mkdir -p ${mapdir}/${label}/star1
fi

if [ ! -e ${mapdir}/${label}/stargenome ]; then
    mkdir -p ${mapdir}/${label}/stargenome
fi

if [ ! -e ${mapdir}/${label}/star2 ]; then
    mkdir -p ${mapdir}/${label}/star2
fi


echo "[$(date)] Start 1st STAR mapping"
# first run
cd ${mapdir}/${label}/star1

STAR --genomeDir ${starrefdir} \
    --readFilesIn ${fqdir}/${fq1} ${fqdir}/${fq2} \
    --readFilesCommand zcat \
    --runThreadN ${thread}

echo "[$(date)] End 1st STAR mapping"
echo "[$(date)] Start STAR build index"
# new STAR index
ln -s ${thefa} ${mapdir}/${label}/stargenome/
samtools faidx ${mapdir}/${label}/stargenome/genome.fa

STAR --runMode genomeGenerate --genomeDir ${mapdir}/${label}/stargenome \
     --genomeFastaFiles ${mapdir}/${label}/stargenome/genome.fa \
     --sjdbFileChrStartEnd ${mapdir}/${label}/star1/SJ.out.tab \
     --sjdbOverhang 75 --runThreadN ${thread}

echo "[$(date)] End STAR build index"
# 2nd alignment
echo "[$(date)] Start 2nd STAR mapping"
cd ${mapdir}/${label}/star2

STAR --genomeDir ${mapdir}/${label}/stargenome \
    --readFilesIn ${fqdir}/${fq1} ${fqdir}/${fq2} \
    --readFilesCommand zcat \
    --runThreadN ${thread}

echo "[$(date)] End 2nd STAR mapping"

# add read groups sort
echo "[$(date)] Start Picard AddOrReplaceReadGroups"

$Picard AddOrReplaceReadGroups \
    I=Aligned.out.sam \
    O=rg_added_sorted.bam \
    SO=coordinate \
    RGID=${label} RGLB=library RGPL=illumina RGPU=xten RGSM=${label}

rm -r Aligned.out.sam
echo "[$(date)] End Picard AddOrReplaceReadGroups"

# setting MAPQ to 0 for unmapped reads
echo "[$(date)] Start Picard CleanSam"

$Picard CleanSam I=rg_added_sorted.bam O=dedupped.bam

echo "[$(date)] End Picard CleanSam"
echo "[$(date)] Start Picard MarkDuplicates"

$Picard MarkDuplicates \
    I=clean.bam\
    O=dedupped.bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    M=output.metrics

echo "[$(date)] End Picard MarkDuplicates"

echo "[$(date)] Start GATK SplitNCigarReads"

gatk SplitNCigarReads \
    -R ${starfa} \
    -I dedupped.bam \
    -O ${mapdir}/${label}.star.bam

echo "[$(date)] End GATK SplitNCigarReads"

rm -fr ${mapdir}/${label}

gzip ${fqdir}/${label}.1.fq
gzip ${fqdir}/${label}.2.fq

####################
echo "[$(date)] ===== ALL FINISHED ====="
####################
