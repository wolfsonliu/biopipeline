#! /bin/bash
#SBATCH --job-name=P_cv
#SBATCH --output=${HOME}/log/%x.%j_%A_%a.%N.array.out
#SBATCH --error=${HOME}/log/%x.%j_%A_%a.%N.array.err
#SBATCH --partition=C032M0128G
#SBATCH --qos=low
#SBATCH --get-user-env
#SBATCH -n 4
#SBATCH --cpu-freq=high
#SBATCH -A hpc
#SBATCH --mail-type=end
#SBATCH --mail-user=
#SBATCH --time=120:00:00

####################
parr=4

         # 0      1       2      3      4           5             6        7        8         9     10      11     12       13
labels=(FANCC_1 FANCC_2 Mock_1 Mock_2 Untreated_1 Untreated_2 PPIB111_1 PPIB111_2 PPIB111_3 PPIB_1 PPIB_2 PPIB_3 Random_1 Random_2)

label=${labels[${SLURM_ARRAY_TASK_ID}]}

####################

basedir=${HOME}/Project
codedir=${basedir}/code
datadir=${basedir}/data
fqdir=${datadir}/clean
mapdir=${datadir}/bam
starrefdir=${datadir}/star
starfa=${starrefdir}/genome.fa
vcfdir=${datadir}/vcf

# Reference
genebed=${datadir}/genegood.bed

refdir=${homedir}/Reference/Homo_sapiens
faref=${refdir}/UCSC/hg38/Sequence/BWAIndex/genome.fa

g1000=${refdir}/1000genome/1000G_phase1.snps.high_confidence.hg38.vcf.gz
hapmap=${refdir}/hapmap/hapmap_3.3.hg38.vcf.gz
omni=${refdir}/omni/1000G_omni2.5.hg38.vcf.gz
dbsnp=${refdir}/dbSNP/dbsnp_146.hg38.vcf.gz # version 146 GRCh38.p7
evs=${refdir}/EVS/EVS.new.vcf.gz
ANNOVARDB=${HOME}/Software/annovar/humandb

####################


if [ ! -e ${vcfdir} ]; then
    mkdir -p ${vcfdir}
fi

####################


outdir=${vcfdir}

# HaplotypeCaller
echo "[$(date)] Start GATK HaplotypeCaller <- "${label}.star.bam
gatk --java-options "-Xmx4g" HaplotypeCaller  \
    -R ${faref} \
    -I ${mapdir}/${label}.star.bam \
    -O ${outdir}/${label}.vcf.gz \
    -L ${genebed} \
    -D ${dbsnp} \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20.0
echo "[$(date)] End GATK HaplotypeCaller -> "${label}.vcf.gz


# hard filter
echo "[$(date)] Start GATK VariantFiltration <- "${label}.rc.vcf.gz
gatk VariantFiltration \
    -R ${faref} \
    -V ${outdir}/${label}.vcf.gz \
    -O ${outdir}/${label}.filtered.vcf.gz \
    -L ${genebed} \
    --cluster-window-size 35 --cluster-size 3 \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP < 20.0 || QUAL < 20" \
    --filter-name "BADSNP"
tabix -f ${outdir}/${label}.filtered.vcf.gz
echo "[$(date)] End GATK VariantFiltration -> "${label}.filtered.vcf.gz

# select SNP
echo "[$(date)] Start GATK SelectVariants <- "${label}.filtered.vcf.gz
gatk SelectVariants \
    -R ${faref} \
    --select-type-to-include SNP \
    -select "FILTER != 'BADSNP'" \
    --discordance ${dbsnp} \
    -V ${outdir}/${label}.filtered.vcf.gz \
    -O ${outdir}/${label}.rmdbsnp.vcf.gz \
    -L ${genebed}
tabix -f ${outdir}/${label}.rmdbsnp.vcf.gz

gatk SelectVariants \
    -R ${faref} \
    --select-type-to-include SNP \
    -select "QUAL > 20" \
    --discordance ${g1000} \
    -V ${outdir}/${label}.rmdbsnp.vcf.gz \
    -O ${outdir}/${label}.rmg1000.vcf.gz \
    -L ${genebed}
tabix -f ${outdir}/${label}.rmg1000.vcf.gz

gatk SelectVariants \
    -R ${faref} \
    --select-type-to-include SNP \
    -select "QUAL > 20" \
    --discordance ${evs} \
    -V ${outdir}/${label}.rmg1000.vcf.gz \
    -O ${outdir}/${label}.snp.vcf.gz \
    -L ${genebed}
tabix -f ${outdir}/${label}.snp.vcf.gz

echo "[$(date)] End GATK SelectVariants -> "${label}.snp.vcf.gz


# annotate vcf with strand and geneid
echo "[$(date)] Start bcftools annotate <- "${label}.snp.vcf.gz
bcftools annotate \
    -Oz -a ${genebed}.gz \
    -c CHROM,FROM,TO,INFO/geneid,-,INFO/strand \
    -h ${datadir}/vcfbed.hdr \
    ${outdir}/${label}.snp.vcf.gz > ${outdir}/${label}.snp.anno.vcf.gz
tabix -f ${outdir}/${label}.snp.anno.vcf.gz
echo "[$(date)] End bcftools annotate -> "${label}.snp.anno.vcf.gz



bcftools filter -Ov -i 'FILTER !~ "BADSNP" && INFO/DP > 50  && QUAL > 20 && ((strand=="+" && REF=="A" && ALT=="G") || (strand=="-" && REF=="T" && ALT=="C"))' ${outdir}/${label}.snp.anno.vcf.gz > ${outdir}/${label}.ag.vcf

# annotate with annovar


table_annovar.pl ${outdir}/${label}.ag.vcf ${ANNOVARDB} \
    --buildver hg38 --remove \
    --outfile ${outdir}/${label}.anv \
    --protocol refGene,exac03,avsnp150 --operation g,f,f \
    --nastring . --vcfinput
####################
echo "[$(date)] ===== ALL FINISHED ====="
####################
