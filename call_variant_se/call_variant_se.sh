#! /bin/bash
function usage {
    echo "Usage: $0 -r reference -f input.fastq -d outputdirectory -l outputfilelabel -t thread" 1>&2
}


while getopts "hr:f:d:l:t:" opt; do
    case $opt in
        h)
            usage
            exit 0
            ;;
        r)
            referencefile=$OPTARG
            ;;
        f)
            # input fastq file
            inputfq=$OPTARG
            ;;
        d)
            # output directory
            outputdir=$OPTARG
            ;;
        l)
            # output name
            outputlabel=$OPTARG
            ;;
        t)
            # thread
            PARR=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

shift $((OPTIND-1))

outdir=${outputdir}
outlabel=${outputlabel}
outprefix=${outdir}/${outlabel}

# read reference file path
source $referencefile


################################################################################


if [ ! -e ${outdir}/log ]; then
    mkdir -p ${outdir}/log
fi


if [ ! -e ${outdir}/bam ]; then
    mkdir -p ${outdir}/bam
fi


if [ ! -e ${outdir}/vcf ]; then
    mkdir -p ${outdir}/vcf
fi

####################
# Quality Control: fastqc
#     input:  fastq
#     output: html
####################

echo "[$(date)] < Start >  FastQC: "${inputfq}
$Fastqc -o ${outdir}/log/ -f fastq ${inputfq} 2> ${outdir}/log/${outlabel}.FastQC.log
if [ ! $? ]; then
    echo "{ Error! }: [$(date)] FastQC"
    cat ${outdir}/log/${outlabel}.FastQC.log
    exit 1
fi
echo "[$(date)] <  End  >  FastQC: result in "${outdir}/log
echo ""

####################
# Mapping: bwa
#     input:  fastq
#     input:  bowtie2index
#     output: sorted.bam
####################

echo "[$(date)] < Start >  bwa: "${inputfq}
$Bwa mem -M ${bwaindexdir}/genome.fa ${inputfq} 2> ${outdir}/log/${outlabel}.bwa.log | \
    $Samtools view -bS - | \
    $Samtools sort -@ ${PARR} - -o ${outdir}/bam/${outlabel}.bam
if [ ! $? ]; then
    echo "{ Error! }: [$(date)] BWA MEM"
    cat ${outdir}/log/${outlabel}.bwa.log
    exit 1
fi
echo "[$(date)] <  End  >  bwa: "${outdir}/bam/${outlabel}.bam
echo ""

####################
# Read Group: picard
#    input:  bam
#    output: sorted.bam
####################

echo "[$(date)] < Start >  Picard AddOrReplaceReadGroups: "${outdir}/bam/${outlabel}.bam
$Picard AddOrReplaceReadGroups \
     I=${outdir}/bam/${outlabel}.bam \
     O=${outdir}/bam/${outlabel}.sorted.bam \
     SORT_ORDER=coordinate \
     RGID=${outlabel} \
     RGLB=bwa \
     RGPL=illumina \
     RGPU=unit1 \
     RGSM=20 2> ${outdir}/log/${outlabel}.picard.addgroup.log
if [ ! $? ]; then
    echo "{ Error! }: [$(date)] Picard AddOrReplaceReadGroups"
    cat ${outdir}/log/${outlabel}.picard.addgroup.log
    exit 1
fi
echo "[$(date)] <  End  >  Picard AddOrReplaceReadGroups: "${outdir}/bam/${outlabel}.sorted.bam
echo ""

####################
# Mark Duplicates: picard
#    input:  sorted.bam
#    output: markdup.sorted.bam
#    output: markdup.txt
####################

echo "[$(date)] < Start >  Picard MarkDuplicates: "${outdir}/bam/${outlabel}.sorted.bam
$Picard MarkDuplicates I=${outdir}/bam/${outlabel}.sorted.bam \
       O=${outdir}/bam/${outlabel}.sorted.markdup.bam \
       M=${outdir}/bam/${outlabel}.sorted.markdup.txt 2> ${outdir}/log/${outlabel}.picard.markduplicates.log
if [ ! $? ]; then
    echo "{ Error! }: [$(date)] Picard MarkDuplicates"
    cat ${outdir}/log/${outlabel}.picard.markduplicates.log
    exit 1
fi
echo "[$(date)] <  End  >  Picard MarkDuplicates: "${outdir}/bam/${outlabel}.sorted.markdup.bam
echo ""

####################
# Statistics of bam file: samtools
#    input:  markdup.sorted.bam
#    output: markdup.sorted.bam.stats
####################

$Samtools stats ${outdir}/bam/${outlabel}.sorted.markdup.bam > ${outdir}/log/${outlabel}.sorted.markdup.bam.stats
if [ ! $? ]; then
    echo "{ Error! }: [$(date)] samtools stats"
    cat ${outdir}/log/${outlabel}.sorted.markdup.bam.stats
    exit 1
fi

####################
# Build bai: samtools
#    input:  markdup.sorted.bam
#    output: markdup.sorted.bai
####################

$Samtools index -b ${outdir}/bam/${outlabel}.sorted.markdup.bam
if [ ! $? ]; then
    echo "{ Error! }: [$(date)] samtools index"
    exit 1
fi

####################
# Call variants: GATK HaplotypeCaller
#    input:  markdup.sorted.bam
#    input:  reference fa
#    output: vcf
####################

echo "[$(date)] < Start >  GATK HaplotypeCaller: "${outdir}/bam/${outlabel}.sorted.markdup.bam
$Gatk --java-options "-Xmx8G" HaplotypeCaller \
     -ERC GVCF \
     -R ${genomefa} \
     --output-mode EMIT_VARIANTS_ONLY \
     -I ${outdir}/bam/${outlabel}.sorted.markdup.bam \
     -O ${outdir}/vcf/${outlabel}.g.vcf.gz 2> ${outdir}/log/${outlabel}.gatk.HaplotypeCaller.log
if [ ! $? ]; then
    echo "{ Error! }: [$(date)] GATK HaplotypeCaller"
    cat ${outdir}/log/${outlabel}.gatk.HaplotypeCaller.log
    exit 1
fi
echo "[$(date)] <  End  >  GATK HaplotypeCaller: "${outdir}/vcf/${outlabel}.g.vcf.gz
echo ""

####################
# Convert gvcf to vcf: GATK GenotypeGVCFs
#    input:  gatk.haplotype.g.vcf
#    output: vcf
####################

echo "[$(date)] < Start >  GATK GenotypeGVCFs: "${outdir}/vcf/${outlabel}.g.vcf.gz
$Gatk --java-options "-Xmx8G" GenotypeGVCFs \
     -R ${genomefa} \
     -V ${outdir}/vcf/${outlabel}.g.vcf.gz \
     -O ${outdir}/vcf/${outlabel}.vcf.gz 2> ${outdir}/log/${outlabel}.gatk.GenotypeGVCFs.log
if [ ! $? ]; then
    echo "{ Error! }: [$(date)] GATK GenotypeGVCFs"
    cat ${outdir}/log/${outlabel}.gatk.GenotypeGVCFs.log
    exit 1
fi
echo "[$(date)] <  End  >  GATK GenotypeGVCFs: "${outdir}/vcf/${outlabel}.vcf.gz
echo ""

####################
# Filter: GATK VariantFiltration
#    input: vcf
#    ouput: vcf
####################

echo "[$(date)] < Start >  GATK VariantFiltration: "${outdir}/vcf/${outlabel}.vcf.gz
$Gatk VariantFiltration \
     -R ${genomefa} \
     -V ${outdir}/vcf/${outlabel}.vcf.gz \
     -O ${outdir}/vcf/${outlabel}.filter0.vcf.gz \
     --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
     --filter-name "SNPfilter" 2> ${outdir}/log/${outlabel}.gatk.VariantFiltration.snp.log
if [ ! $? ]; then
    echo "{ Error! }: [$(date)] GATK VariantFiltration"
    exit 1
fi
echo "[$(date)] <  End  >  GATK VariantFiltration: "${outdir}/bam/${outlabel}.filter.vcf.gz
echo ""

zcat ${outdir}/vcf/${outlabel}.filter0.vcf.gz | grep -v "SNPfilter" > ${outdir}/vcf/${outlabel}.filter.vcf

rm -f ${outdir}/vcf/${outlabel}.filter0.vcf.gz
rm -f ${outdir}/vcf/${outlabel}.filter0.vcf.gz.tbi
####################
# Variants Annotation: annovar
#    input:  vcf.gz
#    output: output dir vcfs
####################


if [ ! -e ${outdir}/humandb ]; then
    echo "[$(date)] < Start >  ANNOVAR annotate_variation.pl: Downloading database from websites"
    # annotate_variation.pl -buildver hg38 --downdb cytoBand ls${outdir}/humandb/ >> ${outdir}/log/${outlabel}.annotate_variation.log
    $Annotate_variation --buildver hg38 --downdb --webfrom annovar refGene ${outdir}/humandb/ >> ${outdir}/log/${outlabel}.annotate_variation.log
    $Annotate_variation --buildver hg38 --downdb --webfrom annovar exac03 ${outdir}/humandb/ >> ${outdir}/log/${outlabel}.annotate_variation.log
    $Annotate_variation --buildver hg38 --downdb --webfrom annovar avsnp147 ${outdir}/humandb/ >> ${outdir}/log/${outlabel}.annotate_variation.log
    $Annotate_variation --buildver hg38 --downdb --webfrom annovar dbnsfp30a ${outdir}/humandb/ >> ${outdir}/log/${outlabel}.annotate_variation.log
    $Annotate_variation --buildver hg38 --downdb --webfrom annovar clinvar_20180603 ${outdir}/humandb/ >> ${outdir}/log/${outlabel}.annotate_variation.log
    echo "[$(date)] <  End  >  ANNOVAR annotate_variation.pl: Downloaded database in " ${outdir}/humandb
    echo ""
fi

echo "[$(date)] < Start >  ANNOVAR table_annovar.pl: "${outdir}/vcf/${outlabel}.filter.vcf
$Table_annovar --buildver hg38 \
              --remove --protocol refGene,clinvar_20180603,exac03,avsnp147,dbnsfp30a \
              --operation g,f,f,f,f -nastring . --polish --vcfinput \
              --out ${outdir}/vcf/${outlabel} \
              ${outdir}/vcf/${outlabel}.filter.vcf ${outdir}/humandb/
echo "[$(date)] <  End  >  ANNOVAR table_annovar.pl: " ${outdir}/vcf/${outlabel}.hg38_multianno.txt
echo ""

cp ${outdir}/vcf/${outlabel}.hg38_multianno.vcf ${outdir}/${outlabel}.anno.vcf

echo "All finished, results in: " ${outdir}
echo "Annotated VCF: " ${outdir}/${outlabel}.anno.vcf
#####################
