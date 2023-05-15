#!/bin/bash

# May 11 2023

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -f <fastq1> -r <fastq2> -a <forward_adapt> -A <reverse_adapt> -p <Prefix> -o <outputdir> -b <blacklist> -m <min_length> -x <index> -n <normalize_method> -s <step>"
	echo ""
	echo " -f	file         	[required] fastq1"
	echo ""
	echo " -r	file         	[required] fastq2"
	echo ""
	echo " -a	string       	[required] forward_adapt"
	echo ""
	echo " -A	string       	[required] reverse_adapt"
	echo ""
	echo " -p	string       	[required] output prefix"
	echo ""
	echo " -o	dir          	[required] output dir"
	echo ""
	echo " -b	file         	[optional] blacklist default:mm10-blacklist.v2.bed"
	echo ""
	echo " -m	number         	[optional] min_length, defult:25"
	echo ""
	echo " -x	path and Prefix [optional] index using for mapping defult:mm10_dm6_bowtie2_index"
	echo ""
	echo " -n	string         	[optional] normalize_method, defult:RPKM;RPKM, CPM, BPM, RPGC, None"
	echo ""
	echo " -s	string         	[optional] which step you want to run: QC,cutadapt,mapping,bigwig,callpeak"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}

min_length=25
index=mm10_dm6_bowtie2_index
normalize_method=RPKM
blacklist=mm10-blacklist.v2.bed

if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "f:r:a:A:p:o:b:m:x:n:s:h" optionName
do
	case $optionName in
		f) fastq1="$OPTARG";;
		r) fastq2="$OPTARG";;
		a) forward_adapt="$OPTARG";;
		A) reverse_adapt="$OPTARG";;
		p) Prefix="$OPTARG";;
		o) outputdir="$OPTARG";;
		b) blacklist="$OPTARG";;
		m) min_length="$OPTARG";;
		x) index="$OPTARG";;
		n) normalize_method="$OPTARG";;
		s) step="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $fastq1 = "" ]]; then
	echo "the $fastq1 file is needed "
	exit 1
elif [[ ! -f $fastq1 ]]; then
	echo "$fastq1:   is not found"
	exit 2
fi

if [[ $fastq2 = "" ]]; then
	echo "the $fastq2 file is needed "
	exit 1
elif [[ ! -f $fastq2 ]]; then
	echo "$fastq2:   is not found"
	exit 2
fi

if [[ $forward_adapt = "" ]]; then
	echo " the $forward_adapt string is needed "
	exit 1
fi

if [[ $reverse_adapt = "" ]]; then
	echo " the $reverse_adapt string is needed "
	exit 1
fi

if [[ $Prefix = "" ]]; then
	echo " the $Prefix string is needed "
	exit 1
fi

if [[ $outputdir = "" ]]; then
	echo "the $outputdir file is needed "
	exit 1
elif [[ ! -d $outputdir ]]; then
#	 echo "$outputdir:   is not found"
#	exit 2
	mkdir -p ${outputdir}
fi

if [[ $blacklist = "" ]]; then
	echo "the $blacklist file is needed "
	exit 1
elif [[ ! -f $blacklist ]]; then
	echo "$blacklist:   is not found"
	exit 2
fi


# Step 1: Trim adapters from paired-end FASTQ files.
# Note: Make sure cutadapt is installed before running this script.

if [[ ! -d ${outputdir}/cutadapt ]]; then
	mkdir -p ${outputdir}/cutadapt/QC
        cutadapt -a ${forward_adapt} -A ${reverse_adapt} -q 15,15 --overlap 1 -m ${min_length} -o ${outputdir}/cutadapt/${Prefix}_1.fq -p ${outputdir}/cutadapt/${Prefix}_2.fq ${fastq1} ${fastq2}

# Step 2: Run FastQC to check the quality of the trimmed reads
# Note: Step 2 is part of a large script with the previous step.
# Variables are defined earlier in the same script.
# run FastQC on the cutadapt output files
	fastqc ${outputdir}/cutadapt/${Prefix}_1.fq ${outputdir}/cutadapt/${Prefix}_2.fq -o ${outputdir}/cutadapt/QC
fi

#--- 3. mapping
if [[ ! -d ${outputdir}/mapping ]]; then
	# create directory for mapping results
	mkdir -p ${outputdir}/mapping
	echo "mapping to genome..."
	# map reads to reference and spike-in genome using Bowtie2
        bowtie2 -q -x ${index} --reorder -1 ${outputdir}/cutadapt/${Prefix}_1.fq -2 ${outputdir}/cutadapt/${Prefix}_2.fq -p 15 -S ${outputdir}/mapping/${Prefix}.sam 2> ${outputdir}/mapping/${Prefix}_mapping.stat
        # Step 6: sort the SAM file and convert it to BAM format
	samtools sort -@ 10 ${outputdir}/mapping/${Prefix}.sam >> ${outputdir}/mapping/${Prefix}.bam
	rm -f ${outputdir}/mapping/${Prefix}.sam
	# Step 6: mark duplicates using sambamba
	sambamba markdup -r -t 10 ${outputdir}/mapping/${Prefix}.bam ${outputdir}/mapping/${Prefix}_rmdup.bam
fi

echo "extracting unique & concordant pairs..."
samtools view -@ 5 -hF 4 ${outputdir}/mapping/${Prefix}_rmdup.bam | grep -v "XS:" | samtools view -@ 5 -bS -o ${outputdir}/mapping/${Prefix}_unique.bam
samtools view -H ${outputdir}/mapping/${Prefix}_rmdup.bam > ${outputdir}/mapping/header.txt
samtools view -@ 5 -hF 4 ${outputdir}/mapping/${Prefix}_unique.bam | grep "YT:Z:CP" | cat ${outputdir}/mapping/header.txt - | samtools view -@ 5 -bS -o ${outputdir}/mapping/${Prefix}_unique_concordant.bam
samtools index ${outputdir}/mapping/${Prefix}_unique_concordant.bam

# Step 7: split the "composite" BAM files into mouse and Drosophila BAM files
samtools view -h ${outputdir}/mapping/${Prefix}_unique_concordant.bam `echo dm6_chr{2L,2R,3L,3R,4,X,Y,M}` -o ${outputdir}/mapping/${Prefix}_unique_concordant_dm6.bam
samtools view -h ${outputdir}/mapping/${Prefix}_unique_concordant.bam `echo chr{{1..19},X,Y,M}` -o ${outputdir}/mapping/${Prefix}_unique_concordant_mm10.bam


# Step 9: Use SAMtools to assign plus and minus strands for uniquely mapped reads
mkdir -p ${outputdir}/mapping/split_strand
samtools view -@ 10 -hf 80 -F 32 ${outputdir}/mapping/${Prefix}_unique_concordant_mm10.bam -o ${outputdir}/mapping/split_strand/${Prefix}_r1_rev.bam
samtools view -@ 10 -hf 96 -F 16 ${outputdir}/mapping/${Prefix}_unique_concordant_mm10.bam -o ${outputdir}/mapping/split_strand/${Prefix}_r1_forw.bam
samtools view -@ 10 -hf 144 -F 32 ${outputdir}/mapping/${Prefix}_unique_concordant_mm10.bam -o ${outputdir}/mapping/split_strand/${Prefix}_r2_rev.bam
samtools view -@ 10 -hf 160 -F 16 ${outputdir}/mapping/${Prefix}_unique_concordant_mm10.bam -o ${outputdir}/mapping/split_strand/${Prefix}_r2_forw.bam

# Step 9: Use SAMtools to assign plus and minus strands for uniquely mapped reads
samtools merge -f -@ 5 ${outputdir}/mapping/split_strand/${Prefix}_plus.bam ${outputdir}/mapping/split_strand/${Prefix}_r1_rev.bam ${outputdir}/mapping/split_strand/${Prefix}_r2_forw.bam
samtools sort -@ 5 -o ${outputdir}/mapping/split_strand/${Prefix}_plus_sorted.bam ${outputdir}/mapping/split_strand/${Prefix}_plus.bam 
samtools index ${outputdir}/mapping/split_strand/${Prefix}_plus_sorted.bam
rm ${outputdir}/mapping/split_strand/${Prefix}_plus.bam

samtools merge -f -@ 5 ${outputdir}/mapping/split_strand/${Prefix}_minus.bam ${outputdir}/mapping/split_strand/${Prefix}_r1_forw.bam ${outputdir}/mapping/split_strand/${Prefix}_r2_rev.bam
samtools sort -@ 5 -o ${outputdir}/mapping/split_strand/${Prefix}_minus_sorted.bam ${outputdir}/mapping/split_strand/${Prefix}_minus.bam
samtools index ${outputdir}/mapping/split_strand/${Prefix}_minus_sorted.bam
rm ${outputdir}/mapping/split_strand/${Prefix}_minus.bam

# Step 10: Generate bigWig files for plus and minus strand reads
# scale_factor=$(echo 1000000/${fly_num} | bc -l) Count reads originating from Drosophila S2 spike-in cells and calculate the calibration factors (alpha=1e6/dm6_count) for reads that mapped to the mouse genome.
if [[ ! -d ${outputdir}/bigwig ]]; then
	mkdir -p ${outputdir}/bigwig
	echo "making bigwig tracks..."
	sort -k1,1 -k2,2n ${blacklist} | bedtools merge -i - > ${outputdir}/bigwig/${Prefix}_blacklist.bed
	bamCoverage --bam ${outputdir}/mapping/split_strand/${Prefix}_plus_sorted.bam -o ${outputdir}/bigwig/${Prefix}_spikein.bw --binSize 10 --normalizeUsing None --scaleFactor ${scale_factor} --effectiveGenomeSize 2652783500 -p 5
        bamCoverage --bam ${outputdir}/mapping/split_strand/${Prefix}_minus_sorted.bam -o ${outputdir}/bigwig/${Prefix}_spikein.bw --binSize 10 --normalizeUsing None --scaleFactor ${scale_factor} --effectiveGenomeSize 2652783500 --effectiveGenomeSize 2652783500 -p 5

fi


# Step 16: Record alignment statistics
if [[ ! -f ${outputdir}/${Prefix}_stastic.csv ]]; then
	echo "calculating mapping stats..."
	pwd
	suffix=$(echo ${fastq1} | awk -F "." '{print $NF}')
	if [[ $suffix == 'gz' ]]; then
		total_reads=$(grep "^@" <(zcat ${fastq1}) | wc -l)
	elif [[ $suffix != 'gz' ]]; then
		total_reads=$(grep "^@" ${fastq1} | wc -l)
	fi

	cutadapt=$(grep "^@" ${outputdir}/cutadapt/${Prefix}_1.fq | wc -l)
	mapping=$(samtools view -F 4 ${outputdir}/mapping/${Prefix}.bam | wc -l)
	rmdup=$(samtools view -F 4 ${outputdir}/mapping/${Prefix}_rmdup.bam | wc -l)
	unique=$(samtools view -c ${outputdir}/mapping/${Prefix}_unique.bam)
	concordant=$(samtools view -c ${outputdir}/mapping/${Prefix}_unique_concordant.bam)
	mm10_concord=$(samtools view -c ${outputdir}/mapping/${Prefix}_unique_concordant_mm10.bam)
	dm6_concord=$(samtools view -c ${outputdir}/mapping/${Prefix}_unique_concordant_dm6.bam)
#--- stastic file
	echo -e "sample_name,total_reads,cutadapt,mapping reads,rm duplicate,unique mapping,concordant pairs,mm10 concordant pairs,dm6 concordant pairs" >> ${outputdir}/${Prefix}_stastic.csv
	echo -e "${Prefix},${total_reads},${cutadapt},${mapping},${rmdup},${unique},${concordant},${mm10_concord},${dm6_concord}" >> ${outputdir}/${Prefix}_stastic.csv
fi


