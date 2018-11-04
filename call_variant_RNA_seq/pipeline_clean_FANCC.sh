#! /bin/bash
#SBATCH --job-name=P_Clean_FANCC
#SBATCH --output=${HOME}/log/%x.%j_%A_%a.%N.out
#SBATCH --error=${HOME}/log/%x.%j_%A_%a.%N.err
#SBATCH --partition=C032M0128G
#SBATCH --qos=low
#SBATCH --get-user-env
#SBATCH -n 1
#SBATCH --cpu-freq=high
#SBATCH -A hpc
#SBATCH --mail-type=end
#SBATCH --mail-user=
#SBATCH --time=120:00:00


# clean contamination
function fq2rmpattern {
    pattern=$1
    file1=$2
    file2=$3

    paste -d '\t' $file1 $file2 | awk -v PATTERN=$pattern 'BEGIN {
        FS="\t";
        havepattern = 0;
    }
    FNR % 4 == 1 {
        a1 = $1;
        a2 = $2;
    }
    FNR % 4 == 2 {
        b1 = $1;
        b2 = $2;
        if ($0 ~ PATTERN) {
            havepattern = 1;
        }
    }
    FNR % 4 == 3 {
        c1 = $1;
        c2 = $2;
    }
    FNR % 4 == 0 {
        d1 = $1;
        d2 = $2;
        if (havepattern == 0) {
            print a1 "\n" b1 "\n" c1 "\n" d1 > "/dev/stdout";
            print a2 "\n" b2 "\n" c2 "\n" d2 > "/dev/stderr";
        }
        havepattern = 0
    }' -
}


function fq2getpattern {
    pattern=$1
    file1=$2
    file2=$3
    paste -d '\t' $file1 $file2 | awk -v PATTERN=$pattern 'BEGIN {
        FS="\t";
        havepattern = 0;
    }
    FNR % 4 == 1 {
        a1 = $1;
        a2 = $2;
    }
    FNR % 4 == 2 {
        b1 = $1;
        b2 = $2;
        if ($0 ~ PATTERN) {
            havepattern = 1;
        }
    }
    FNR % 4 == 3 {
        c1 = $1;
        c2 = $2;
    }
    FNR % 4 == 0 {
        d1 = $1;
        d2 = $2;
        if (havepattern == 1) {
            print a1 "\n" b1 "\n" c1 "\n" d1 > "/dev/stdout";
            print a2 "\n" b2 "\n" c2 "\n" d2 > "/dev/stderr";
        }
        havepattern = 0
    }' -
}


####################


Pcontrol="TCTCAGTCCAATGTATGGTCCGAGCACAAGCTCTAATCAAAGTCCGCGGGTGTAGACCGGTTGCCATAGGA"
Qcontrol=$(echo ${Pcontrol} | rev | tr "ATGC" "TACG")
PguideHead="GCCTCCCATCACGGGGGCCGT"
QguideHead="ACGGCCCCCGTGATGGGAGGC"
PguideTail="TCAAAGGGACCTCCGCAGTTTT{1,4}"
QguideTail="A{1,4}AAACTGCGGAGGTCCCTTTGA"

basedir=${HOME}/Project
datadir=${basedir}/data
rawdatadir=${datadir}/rawfq
cleandir=${datadir}/clean
qcdir=${datadir}/qc
guidedir=${datadir}/guide
tmp1dir=${datadir}/tmp1
tmp2dir=${datadir}/tmp2

Fastqc=${HOME}/Software/FastQC/fastqc
####################

if [ ! -e ${cleandir} ]; then
    mkdir -p ${cleandir}
fi

if [ ! -e ${qcdir} ]; then
    mkdir -p ${qcdir}
fi

if [ ! -e ${guidedir} ]; then
    mkdir -p ${guidedir}
fi

if [ ! -e ${tmp1dir} ]; then
    mkdir -p ${tmp1dir}
fi


if [ ! -e ${tmp2dir} ]; then
    mkdir -p ${tmp2dir}
fi

filenames=(FAN_1_S2_L008 FAN_2_S3_L008 M_2_S3_L007 M_3_S1_L008 UN_1_S1_L007 UN_2_S2_L007)
filename=${filenames[${SLURM_ARRAY_TASK_ID}]}
labels=(FANCC_1 FANCC_2 Mock_1 Mock_2 Untreated_1 Untreated_2)

label=${labels[${SLURM_ARRAY_TASK_ID}]}


####################
# Unzip fq
echo "["$(date)"] Unzipping file..."
echo "    "${filename}_R1_001.fastq.gz
echo "    "${filename}_R2_001.fastq.gz

gunzip ${rawdatadir}/${filename}_R1_001.fastq.gz
gunzip ${rawdatadir}/${filename}_R2_001.fastq.gz

echo "["$(date)"] File unzipped..."
echo "    "${filename}_R1_001.fastq
echo "    "${filename}_R2_001.fastq

# QC
$Fastqc -o ${qcdir} \
    ${rawdatadir}/${filename}_R1_001.fastq
$Fastqc -o ${qcdir} \
    ${rawdatadir}/${filename}_R2_001.fastq


####################
# PguideHead
echo "["$(date)"] Removing PguideHead..."
fq2rmpattern ${PguideHead} \
    ${rawdatadir}/${filename}_R1_001.fastq \
    ${rawdatadir}/${filename}_R2_001.fastq \
    1> ${tmp1dir}/${label}.rmPH.1.fq \
    2> ${tmp1dir}/${label}.rmPH.2.fq
echo "["$(date)"] Remove PguideHead finished"
echo "    "${label}.rmPH.1.fq
echo "    "${label}.rmPH.2.fq

# QguideHead
echo "["$(date)"] Removing QguideHead..."
fq2rmpattern ${QguideHead} \
    ${tmp1dir}/${label}.rmPH.1.fq \
    ${tmp1dir}/${label}.rmPH.2.fq \
    1> ${tmp2dir}/${label}.rmPHQH.1.fq \
    2> ${tmp2dir}/${label}.rmPHQH.2.fq
echo "["$(date)"] QguideHead finished"
echo "    "${label}.rmPHQH.1.fq
echo "    "${label}.rmPHQH.2.fq

rm -I ${tmp1dir}/${label}.rmPH.1.fq
rm -I ${tmp1dir}/${label}.rmPH.2.fq

# PguideTail
echo "["$(date)"] Removing PguideTail..."
fq2rmpattern ${PguideTail} \
    ${tmp2dir}/${label}.rmPHQH.1.fq \
    ${tmp2dir}/${label}.rmPHQH.2.fq \
    1> ${tmp1dir}/${label}.rmPHQHPT.1.fq \
    2> ${tmp1dir}/${label}.rmPHQHPT.2.fq
echo "["$(date)"] Remove PguideTail finished"
echo "    "${label}.rmPHQHPT.1.fq
echo "    "${label}.rmPHQHPT.2.fq

rm -I ${tmp2dir}/${label}.rmPHQH.1.fq
rm -I ${tmp2dir}/${label}.rmPHQH.2.fq

# QguideTail
echo "["$(date)"] Removing QguideTail..."
fq2rmpattern ${QguideTail} \
    ${tmp1dir}/${label}.rmPHQHPT.1.fq \
    ${tmp1dir}/${label}.rmPHQHPT.2.fq \
    1> ${tmp2dir}/${label}.rmPHQHPTQT.1.fq \
    2> ${tmp2dir}/${label}.rmPHQHPTQT.2.fq
echo "["$(date)"] Remove PguideTail finished"
echo "    "${label}.rmPHQHPTQT.1.fq
echo "    "${label}.rmPHQHPTQT.2.fq

rm -I ${tmp1dir}/${label}.rmPHQHPT.1.fq
rm -I ${tmp1dir}/${label}.rmPHQHPT.2.fq


# Pcontrol
echo "["$(date)"] Removing Pcontrol..."
fq2rmpattern ${Pcontrol} \
    ${tmp2dir}/${label}.rmPHQHPTQT.1.fq \
    ${tmp2dir}/${label}.rmPHQHPTQT.2.fq \
    1> ${tmp1dir}/${label}.rmPHQHPTQTPC.1.fq \
    2> ${tmp1dir}/${label}.rmPHQHPTQTPC.2.fq
echo "["$(date)"] Remove Pcontrol finished"
echo "    "${label}.rmPHQHPTQTPC.1.fq
echo "    "${label}.rmPHQHPTQTPC.2.fq

rm -I ${tmp2dir}/${label}.rmPHQHPTQT.1.fq
rm -I ${tmp2dir}/${label}.rmPHQHPTQT.2.fq


# Qcontrol
echo "["$(date)"] Removing Qcontrol..."
fq2rmpattern ${Qcontrol} \
    ${tmp1dir}/${label}.rmPHQHPTQTPC.1.fq \
    ${tmp1dir}/${label}.rmPHQHPTQTPC.2.fq \
    1> ${tmp2dir}/${label}.rmPHQHPTQTPCQC.1.fq \
    2> ${tmp2dir}/${label}.rmPHQHPTQTPCQC.2.fq
echo "["$(date)"] Remove Qcontrol finished"
echo "    "${label}.rmPHQHPTQTPCQC.1.fq
echo "    "${label}.rmPHQHPTQTPCQC.2.fq


rm -I ${tmp1dir}/${label}.rmPHQHPTQTPC.1.fq
rm -I ${tmp1dir}/${label}.rmPHQHPTQTPC.2.fq


####################
# cutadapt

forward=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
backward=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
echo "["$(date)"] cutadapt..."
cutadapt -a ${forward} -A ${backward} \
    -u 6 \
    -q 10,20 --minimum-length 90 \
    -o ${cleandir}/${label}.1.fq \
    -p ${cleandir}/${label}.2.fq \
    ${tmp2dir}/${label}.rmPHQHPTQTPCQC.1.fq ${tmp2dir}/${label}.rmPHQHPTQTPCQC.2.fq

rm -I ${tmp2dir}/${label}.rmPHQHPTQTPCQC.1.fq
rm -I ${tmp2dir}/${label}.rmPHQHPTQTPCQC.2.fq
echo "["$(date)"] cutadapt finished."
# QC
$Fastqc -o ${qcdir} \
    ${cleandir}/${label}.1.fq
$Fastqc -o ${qcdir} \
    ${cleandir}/${label}.2.fq

echo "["$(date)"] zipping..."
gzip ${cleandir}/${label}.1.fq
gzip ${cleandir}/${label}.2.fq
echo "["$(date)"] zip finished."
####################
echo "[$(date)] ===== ALL FINISHED ====="
####################
