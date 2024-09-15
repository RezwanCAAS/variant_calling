##########################################
#                                        #
#    fastp for cleaning the raw reads    #
#                                        #
##########################################

#!/bin/bash
#SBATCH --array=1-64
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=8
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --job-name=DraF1_64_filtering
#SBATCH --output=DraF1_64_filtering.%j.out
#SBATCH --mail-user=tariqr@kaust.edu.sa

module load fastp/0.23.2

SEQDIR='~/test/'
VAR=($(ls ${SEQDIR}))
OUTDIR='~/mapped_reads'

N=${SLURM_ARRAY_TASK_ID}

V=${VAR[${N}]}
FASTQ=($(ls ${SEQDIR}${V}/*fastq.gz))
OUTFILE=${OUTDIR}${V}

fastp --in1 ${FASTQ[0]} --in2 ${FASTQ[1]} --out1 ${OUTDIR}${V}_R1_trim.fastq.gz --out2 ${OUTDIR}${V}_R2_trim.fastq.gz -h ${OUTDIR}${V}_trimmed.html -w 30


########################################
#                                      # 
#   BWA for mapping the clean reads    #
#                                      #
########################################


#!/bin/bash
#SBATCH --array=1-64
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=40
#SBATCH --time=300:00:00
#SBATCH --mem-per-cpu=40G
#SBATCH --job-name=DraF1_64_mapping
#SBATCH --output=DraF1_64_mapping.%A_%a.out
#SBATCH --mail-user=tariqr@kaust.edu.sa

module load bwa/0.7.17/gnu-6.4.0 samtools/1.8

SEQDIR='~/draf_wgs/'
VAR=($(ls ${SEQDIR}))
REF='~/genome_data/PitayaGenomic.fa'
OUTDIR='~/mapping_output/'

N=${SLURM_ARRAY_TASK_ID}

V=${VAR[${N}-1]}
FASTQ=($(ls ${SEQDIR}${V}/*trim.fastq.gz))
OUTFILE=${OUTDIR}${V}_Guanhuabai.bam

bwa mem -t 32 ${REF} <(zcat ${FASTQ[0]}) <(zcat ${FASTQ[1]}) | samtools sort -o ${OUTFILE}

########################################
#                                      # 
#      picard for markduplication      #
#           the mapped reads           #
#                                      #
########################################
# After amarkduplication, need to collect the percentage of the duplications. So check the supplementary file

#!/bin/bash
#
#SBATCH --array=1-64
#SBATCH --output=array_markdup.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --mail-user=tariqr@kaust.edu.sa
#SBATCH --mem-per-cpu=40G

module load picard/2.20.4

bam_list=(~/mapped/*.bam)

outputdir=(~/marked/)

bam=${bam_list[$SLURM_ARRAY_TASK_ID-1]}

base=$(basename $bam)

output=${outputdir}${base}

metrics=${outputdir}${base}

picard -Xmx24G MarkDuplicates I=$bam O=${output}_marked.bam  M=${metrics}_marked_metrics.txt ASSUME_SORTED=true COMPRESSION_LEVEL=9 VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE



#################################################
#                                               #
#           GATK for haplotypecalling           #
#                                               #
#                                               #
#################################################


# 1- use this command to add or replace the readgroups. Sometimes bam files are without the readgroups and need to add before haplotypecalling in gatk

#!/bin/bash
#
#SBATCH --array=1-64
#SBATCH --output=GATK_readgroup_1.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mail-user=tariqr@kaust.edu.sa
#SBATCH --mem-per-cpu=30G


module load gatk/4.1.2.0

bam_dir="~/dragon-fruit/markduplicate/"
echo "${bam_dir[@]}"

out_dir="~/dragon-fruit/add_readgroup/"

bam_files=(${bam_dir}/*.bam)
echo "$bam_files"
bam_file=${bam_files[$SLURM_ARRAY_TASK_ID-1]}
base=$(basename $bam_file)
output=${out_dir}${base}

gatk AddOrReplaceReadGroups -I ${base} -O ${output}_RG.bam --RGID readgroup1 --RGLB library1 --RGPL ILLUMINA --RGPU flowcell1_lane1 --RGSM sample1



#2- samtools command to index the addgroup.bam files before haplotypecaller argument in gatk
#!/bin/bash
#
#SBATCH --array=1-64
#SBATCH --output=haplotype_index.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mail-user=tariqr@kaust.edu.sa
#SBATCH --mem-per-cpu=30G

module load samtools/1.8

bam_list=(~/dragon-fruit/add_readgroup/*.bam)

out_dir="~/dragon-fruit/add_readgroup/"

bam_file=${bam_list[$SLURM_ARRAY_TASK_ID-1]}

base=$(basename $bam_file)

output=${out_dir}${base}

samtools index ${bam_file}

#3- haplotypecalling in gatk. put the RG_bam and respective index files in the same directory so that command can go through the index files easily

#!/bin/bash
#
#SBATCH --array=1-64
#SBATCH --output=GATK_variant.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH --mail-user=tariqr@kaust.edu.sa
#SBATCH --mem-per-cpu=30G

module load gatk/4.1.2.0

bam_dir=(~/add_readgroup/*.bam)

out_dir="`~/variant_calling/"

bam_files=${bam_dir[$SLURM_ARRAY_TASK_ID-1]}

base=$(basename $bam_files)

output=${out_dir}${base}

gatk --java-options "-Xmx32g" HaplotypeCaller --native-pair-hmm-threads 20 -R /ibex/scratch/projects/c2141/tariqr/dragonfruits_sequencing/pitaya_genome/PitayaGenomic.fa -I ${bam_files} -O ${output}.g.vcf -ERC GVCF


################################
#                              #
#          CombineGVCF         #
#                              #
###############################

#after haplotyping, the output .g.vcf file will be combined using CombineGVCF tool in gatk


#!/bin/bash
#
#SBATCH --job-name=combine_files
#SBATCH --output=combine_files.%j.out
#SBATCH --err=combine_files.%j.err
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=110:00:00
#SBATCH --mem=800G

module load gatk/4.1.2.0

ref_dir=(~/pitaya_genome/PitayaGenomic.fa)


gatk --java-options "-Xmx800g" CombineGVCFs -R ${ref_dir} --variant vcfs.list -O joint_files.g.vcf.gz


################################
#                              #
#  reheader with bcftools      #
#                              #
###############################

#!/bin/bash
#
#SBATCH --array=1-64
#SBATCH --output=bcftools_header.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH --mail-user=tariqr@kaust.edu.sa
#SBATCH --mem=35G

module load bcftools/1.10.2

bam_dir=(~/variant_calling/*.g.vcf)

out_dir="~/reheading/"

bam_files=${bam_dir[$SLURM_ARRAY_TASK_ID-1]}

base=$(basename $bam_files)

output=${out_dir}${base}

bcftools reheader --samples <(echo "${base}") ${bam_files} > ${output}

################################
#                              #
#  GenotypeGVCF using gatk     #
#                              #
###############################

#genotyped the vcf after combining the vcf files
#!/bin/bash
#
#SBATCH --job-name=genotyped_files
#SBATCH --output=genotyped_files.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=120:00:00
#SBATCH --mem=800G

module load gatk/4.1.2.0

ref_dir=(~/pitaya_genome/PitayaGenomic.fa)

gatk --java-options "-Xmx800g" GenotypeGVCFs -R ${ref_dir} -V combine_variants.g.vcf.gz -O genotyped.vcf.gz





#to play with vcf file to filter the vcf file
 
#gatk for hard filtering the vcf file
https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering

#minor allel frequency
MAF 10%


#to check the number of variants in a unfiltered file

#!/bin/bash
#
#SBATCH --job-name=variants_summary
#SBATCH --output=variants_summary.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=500G

module load bcftools/1.10.2

input_file="genotyped_edit.vcf.gz"
output_file="variants_summary.txt"

bcftools view -H $input_file | wc -l > $output_file

bcftools view -H ~/genotyped.vcf.gz | wc -l > ~/reheading/summary_variants.txt

#to find the minor allele frequency in unfiltered file

#if need to edit the name of the header as per the requirement of the analysis

#!/bin/bash
#
#SBATCH --job-name=reheader_name_edit
#SBATCH --output=reheader_name_edit.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=200G


input_file="genotyped.vcf.gz"
output_file="genotyped_edit.vcf.gz"

#for zipped file
zcat $input_file | sed -e 's/_marked.bam-RG.bam.g.vcf//g' > $output_file
zcat $input_file | sed 's/_marked.bam-RG.bam.g.vcf//' | gzip > $output_file

#for unzip file
sed 's/_marked.bam-RG.bam.g.vcf//' $input_file > $output_file

#to calculate the depth per individual

#!/bin/bash
#
#SBATCH --job-name=depth_individual
#SBATCH --output=depth_individual.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=500G

module load vcftools/0.1.17

input_file="genotyped_edit.vcf.gz"
output_file="depth_individual.txt"

vcftools --gzvcf $input_file --depth --out $output_file



#to calculate the depth per site

#!/bin/bash
#
#SBATCH --job-name=depth_per_site
#SBATCH --output=depth_per_site.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=500G

module load vcftools/0.1.17

input_file="genotyped_edit.vcf.gz"
output_file="depth_per_site.txt"

vcftools --gzvcf $input_file --site-mean-depth --out $output_file


#to calculate site quality

#!/bin/bash
#
#SBATCH --job-name=site_quality
#SBATCH --output=site_quality.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=500G

module load vcftools/0.1.17

input_file="genotyped_edit.vcf.gz"
output_file="site_quality.txt"

vcftools --gzvcf $input_file --site-quality --out $output_file


#to calculate allel frequency

#!/bin/bash
#
#SBATCH --job-name=allel_frequency_2allel
#SBATCH --output=allel_frequency_2allel.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=500G

module load vcftools/0.1.17

input_file="genotyped_edit.vcf.gz"
output_file="allel_frequency_2allel"

vcftools --gzvcf $input_file --freq2 --out $output_file --max-alleles 2


#to calculate allel frequency

#!/bin/bash
#
#SBATCH --job-name=allel_frequency
#SBATCH --output=allel_frequency.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=500G

module load vcftools/0.1.17

input_file="genotyped_edit.vcf.gz"
output_file="allel_frequency"

vcftools --gzvcf $input_file --freq2 --out $output_file


#calculate missing data per individual

#!/bin/bash
#
#SBATCH --job-name=missing_data
#SBATCH --output=missing_data.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=500G

module load vcftools/0.1.17

input_file="genotyped_edit.vcf.gz"
output_file="missing_data.txt"

vcftools --gzvcf $input_file --missing-indv --out $output_file

#calculate the heterozygosity and inbreeding coefficnet per individual 

#!/bin/bash
#
#SBATCH --job-name=individual_heterozygosity
#SBATCH --output=individual_heterozygosity.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=50G

module load vcftools/0.1.17

input_file="genotyped_edit.vcf.gz"
output_file="individual_heterozygosity"

vcftools --gzvcf $input_file --het --out $output_file


#filtering the vcf file with designed parameters
#!/bin/bash
#
#SBATCH --job-name=genotypes_filtering
#SBATCH --output=genotypes_filtering.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=40:00:00
#SBATCH --mem=50G

module load vcftools/0.1.17

input_file="genotyped_edit.vcf.gz"
output_file="genotyped_filtering.vcf.gz"

MAF=0.05
MISS=0.2
QUAL=30
Min_depth=10
Max_depth=40

vcftools --gzvcf $input_file \
         --remove-indels \
         --maf $MAF \
         --max-missing $MISS \
         --minQ $QUAL \
         --min-meanDP $Min_depth \
         --max-meanDP $Max_depth \
         --minDP $Min_depth \
         --maxDP $Max_depth \
         --recode --stdout \
         | gzip -c > $output_file





#calculate the number of variants in vcf after filtering
#!/bin/bash
#
#SBATCH --job-name=filteredvariants_summary
#SBATCH --output=filteredvariants_summary.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=500G

module load bcftools/1.10.2

input_file="genotyped_filtering.vcf.gz"
output_file="filteredvariants_summary.txt"

bcftools view -H $input_file | wc -l > $output_file


#index the filtered vcf file for subset
#!/bin/bash
#
#SBATCH --job-name=indexing_filtered
#SBATCH --output=indexing_filtered.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=40:00:00
#SBATCH --mem=400G

module load vcftools/0.1.17

vcf-sort genotyped_filteringfile2.vcf.gz | bgzip -c > genotyped_filteringfile2_sorted.vcf.gz

tabix -p vcf genotyped_filteringfile2_sorted.vcf.gz


#subset the filtered
#!/bin/bash
#
#SBATCH --job-name=subset_filtering_filtered
#SBATCH --output=subset_filtering.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=40:00:00
#SBATCH --mem=400G

module load gatk/4.1.2.0


input_file="genotyped_filteringfile2_sorted.vcf.gz"
output_file="subset_filtered.vcf.gz"

gatk SelectVariants -V $input_file --select-random-fraction 0.05 -O $output_file


#making hapmap file from subset vcf file

#!/bin/bash
#
#SBATCH --job-name=hapmap_file
#SBATCH --output=hapmap_file.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=40:00:00
#SBATCH --mem=400G

module load plink/1.90b6.24

input_file="subset_filtered.vcf.gz"
output_file="genotyped_hapmap"


plink --vcf $input_file --recode  --const-fid --allow-extra-chr --out $output_file




#convert the vcf file into hapmap file

#!/bin/bash
#
#SBATCH --job-name=hapmap_file
#SBATCH --output=hapmap_file.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=200G

module load vcftools/0.1.17

input_file="subset_filtered.vcf.gz"
output_file="genotyped_hapmap"

vcftools --gzvcf $input_file --plink --out $output_file 


#!/bin/bash
#
#SBATCH --job-name=hapmap_file
#SBATCH --output=hapmap_file.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=200G

module load plink/1.90b6.24

plink --file genotyped_hapmap.map --make-bed --out genotyped_hapmap

plink --file genotyped_hapmap --make-bed --out genotyped_hapmap_binary


#calculate the kinship matrix

#!/bin/bash
#
#SBATCH --job-name=kinship_matrixupdated
#SBATCH --output=kinship_matrixupdated.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=200G


module load vcftools/0.1.17

vcftools --gzvcf subset_filtered.vcf.gz --relatedness2 --out kinship_matrix


#calulcate the genetic distance for phylogenetic tree

#!/bin/bash
#
#SBATCH --job-name=genetic_distance
#SBATCH --output=genetic_distance.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=200G

module load plink/1.90b6.24

plink --file genotyped_hapmap --distance square --out genetic_distance


#calulcate the K values for population structure

#!/bin/bash
#
#SBATCH --job-name=structure_kvalues
#SBATCH --output=structure_kvalues.%j.out
#SBATCH --partition=batch
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=200G

module load admixture/1.3.0

for k in 1 2 3 4 5

do
  admixture --cv genotyped_hapmap.bed $k | tee log${K}.out
done
