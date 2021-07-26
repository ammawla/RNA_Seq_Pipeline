#!/bin/sh

# Args

cores=$1
mem=$2
part=$3
starCores=$4
bamCores=$5
genome=$6
length=$7
type=$8
rawDir=$9
outDir=$10
f=$11
adapter=$12
fastqcAdapters=$13
fastqcContaminants=$14


# Program Paths
qcPath="/share/apps/fastqc-0.11.4/fastqc "
trimPath="/share/apps/Trimmomatic-0.36/trimmomatic-0.36.jar"
starPath="/home/amawla/Data/bin/STAR-2.5/Linux_x86_64_static/STAR"
quantPath="/home/amawla/Data/bin/subread-1.6.3-Linux-x86_64/bin/featureCounts"


if [ $genome = "mm10" ]
then
	GTF="/home/amawla/Data/genomes/Mus_musculus/M8/Annotations/gencode.vM8.primary_assembly.annotation.gtf"
	fasta="/home/amawla/Data/genomes/Mus_musculus/M8/FASTA/GRCm38.primary_assembly.genome.fa"
	if [ $length = "50" ]
	then
		starIndex="/home/amawla/Data/genomes/Mus_musculus/M8/IvsM_STAR_Index_M8/"
	elif [ $length = "100" ]
	then
		starIndex="/home/amawla/Data/genomes/Mus_musculus/M8/STAR_M8_100_Index"
	else
        	echo "Enter a valid read length. Either '50' or '100'"
        	exit 2
	fi
elif [ $genome = "hg19" ]
then
	GTF="/home/amawla/Data/genomes/Homo_sapiens/Gencode/24/gencode.v24.annotation.gtf"
	fasta="/home/amawla/Data/genomes/Homo_sapiens/Gencode/24/FASTA/GRCh38.primary_assembly.genome.fa"
	if [ $length = "50" ]
	then
		starIndex="/home/amawla/Data/genomes/Homo_sapiens/Gencode/24/STAR_Index_24"
	elif [ $length = "100" ]
	then
		starIndex="/home/amawla/Data/genomes/Homo_sapiens/Gencode/24/STAR-100"
	else
        	echo "Enter a valid read length. Either '50' or '100'"
        	exit 2
	fi
else
	echo "Enter an available genome. Either 'mm10' or 'hg19'"
	exit 1
fi

fastQC="$outDir/FastQC"
fastQCBefore="$fastQC/Original"
fastQCAfter="$fastQC/Trimmed"
trimOut="$outDir/Fastq_Trimmed"
starOut="$outDir/BAM"
quantOut="$outDir/featureCounts"
quantGeneOut="$quantOut/gene_level"
quantExonOut="$quantOut/exon_level"
genomeOut="$starOut/Genome_BAM"
transcriptomeOut="$starOut/Transcriptome_BAM"
wigOut="$outDir/Wiggle"
starMetricsOut="$starOut/Metrics"

mkdir -p $fastQC
mkdir -p $fastQCBefore
mkdir -p $fastQCAfter
mkdir -p $trimOut
mkdir -p $starOut
mkdir -p $quantOut
mkdir -p $quantGeneOut
mkdir -p $quantExonOut
mkdir -p $genomeOut
mkdir -p $transcriptomeOut
mkdir -p $wigOut
mkdir -p $starMetricsOut

fastQCScript="$qcPath -t $cores -f fastq -a $fastqcAdapters --contaminants $fastqcContaminants -o "

calc(){ awk "BEGIN { print "$*" }"; }
# Fix later to calculate
minLength=36
trimScript="ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$minLength "

alignScript="$starPath --runMode alignReads -runThreadN $starCores --sjdbOverhang $length --outBAMsortingThreadN $bamCores --genomeDir $starIndex --outWigType wiggle --outWigStrand Unstranded --outSAMattributes NH HI AS NM MD --outFilterMatchNmin 13 --outFilterMismatchNoverReadLmax 0.03 --seedSearchStartLmax 15 --quantMode TranscriptomeSAM GeneCounts --outFilterType BySJout --sjdbInsertSave All --quantTranscriptomeBan IndelSoftclipSingleend --outSAMtype BAM SortedByCoordinate "

quantScript="$quantPath -T $cores -a $GTF -t -G $fasta "

# QC on Original
if [ $type = "single" ]
then
	if [ ! -f "$fastQCBefore/${f}_1_fastqc.html" ];
	then
        	srun -c $cores --partition $part --mem $mem --time 1000 $fastQCScript $fastQCBefore $rawDir/${f}_1.fastq > $fastQCBefore/$f.out 2>&1
	fi
elif [ $type = "paired" ]
then
	if [ ! -f "$fastQCBefore/${f}_1_fastqc.html" ] &&  [ ! -f "$fastQCBefore/${f}_2_fastqc.html" ];
	then
		srun -c $cores --partition $part --mem $mem --time 1000 $fastQCScript $fastQCBefore $rawDir/${f}_1.fastq $rawDir/${f}_2.fastq > $fastQCBefore/$f.out 2>&1
	fi
else
        echo "Enter a valid library type. Either 'single' or 'paired'"
        exit 3
fi


# Trim
if [ $type = "single" ]
then
	if [ ! -f "$trimOut/${f}_1.fastq" ]
	then
	        srun -c $cores --partition $part --mem $mem --time 1000 java -jar $trimPath SE -threads $cores $rawDir/${f}_1.fastq $trimOut/${f}_1.fastq $trimScript > $trimOut/$f.out 2>&1
	fi
elif [ $type = "paired" ]
then
	if [ ! -f "$trimOut/${f}_1P.fastq" ] && [ ! -f "$trimOut/${f}_2P.fastq" ]
	then
		srun -c $cores --partition $part --mem $mem --time 1000 java -jar $trimPath PE  -threads $cores $rawDir/${f}_1.fastq $rawDir/${f}_2.fastq -baseout $trimOut/${f}.fastq $trimScript > $trimOut/$f.out 2>&1
	fi
else
        echo "Enter a valid library type. Either 'single' or 'paired'"
        exit 3
fi

# QC on Trimmed
if [ $type = "single" ]
then
	if [ ! -f "$fastQCAfter/${f}_1_fastqc.html" ]
	then
        	srun -c $cores --partition $part --mem $mem --time 1000 $fastQCScript $fastQCAfter $trimOut/${f}_1.fastq > $fastQCAfter/$f.out 2>&1
	fi
elif [ $type = "paired" ]
then
        if [ ! -f "$fastQCAfter/${f}_1P_fastqc.html" ] &&  [ ! -f "$fastQCAfter/${f}_2P_fastqc.html" ]
	then
	        srun -c $cores --partition $part --mem $mem --time 1000 $fastQCScript $fastQCAfter $trimOut/${f}_1P.fastq $trimOut/${f}_2P.fastq > $fastQCAfter/$f.out 2>&1
	fi
else
        echo "Enter a valid library type. Either 'single' or 'paired'"
        exit 3
fi



# Align
if [ $type = "single" ]
then
	if [ ! -f "$starMetricsOut/${f}Log.final.out" ]
	then
		srun -c $cores --partition $part --mem $mem --time 1000 $alignScript --outFileNamePrefix $starOut/${f%.fastq} --readFilesIn $trimOut/${f}_1.fastq > $starOut/$f.out 2>&1
	fi
elif [ $type = "paired" ]
then
	if [ ! -f "$starMetricsOut/${f}Log.final.out" ]
	then
		srun -c $cores --partition $part --mem $mem --time 1000 $alignScript --outFileNamePrefix $starOut/${f%.fastq} --readFilesIn $trimOut/${f}_1P.fastq $trimOut/${f}_2P.fastq > $starOut/$f.out 2>&1
	fi
else
	echo "Enter a valid library type. Either 'single' or 'paired'"
	exit 3
fi


# Organize bam results
[ -f "$starOut/${f}Aligned.sortedByCoord.out.bam" ] && mv $starOut/${f}Aligned.sortedByCoord.out.bam $genomeOut
[ -f "$starOut/${f}Aligned.toTranscriptome.out.bam" ] && mv $starOut/${f}Aligned.toTranscriptome.out.bam $transcriptomeOut
[ -f "$starOut/${f}Aligned.toTranscriptome.out.bam" ] && mv $starOut/${f}Aligned.toTranscriptome.out.bam $transcriptomeOut
[ -f "$starOut/${f}Log.final.out" ] && mv $starOut/${f}Log.final.out  $starMetricsOut
[ -f "$starOut/${f}Signal.Unique.str1.out.wig" ] && mv $starOut/${f}Signal.Unique.str1.out.wig $wigOut
[ -f "$starOut/${f}_STARtmp" ] && rm -r $starOut/${f}_STARtmp
[ -f "$starOut/${f}Log.out" ] && rm -r $starOut/${f}Log.out
[ -f "$starOut/${f}Log.progress.out" ] && rm -r $starOut/${f}Log.progress.out
[ -f "$starOut/${f}.out" ] && rm -r $starOut/${f}.out
[ -f "$starOut/${f}ReadsPerGene.out.tab" ] && rm -r $starOut/${f}ReadsPerGene.out.tab
[ -f "$starOut/${f}Signal.UniqueMultiple.str1.out.wig" ] && rm -r $starOut/${f}Signal.UniqueMultiple.str1.out.wig
[ -f "$starOut/${f}SJ.out.tab" ] && rm -r $starOut/${f}SJ.out.tab

# Quantify gene level
if [ ! -f "$quantGeneOut/$f.featureCounts.txt" ]
then
	if [ $type = "single" ]
	then
		srun -c $cores --partition $part --mem $mem --time 1000 $quantScript -t exon -g gene_name -o $quantGeneOut/$f.featureCounts.txt $genomeOut/${f}Aligned.sortedByCoord.out.bam  > $quantGeneOut/$f.out 2>&1
	elif [ $type = "paired" ]
	then
		srun -c $cores --partition $part --mem $mem --time 1000 $quantScript -p -B -t exon -g gene_name -o $quantGeneOut/$f.featureCounts.txt $genomeOut/${f}Aligned.sortedByCoord.out.bam  > $quantGeneOut/$f.out 2>&1
	else
        	echo "Enter a valid library type. Either 'single' or 'paired'"
        	exit 3
	fi
fi

# Quantify exon level
if [ ! -f "$quantExonOut/$f.featureCounts.txt" ]
then
	if [ $type = "single" ]
	then
        	srun -c $cores --partition $part --mem $mem --time 1000 $quantScript -t exon -g exon_id -o $quantExonOut/$f.featureCounts.txt $genomeOut/${f}Aligned.sortedByCoord.out.bam  > $quantExonOut/$f.out 2>&1
	elif [ $type = "paired" ]
	then
        	srun -c $cores --partition $part --mem $mem --time 1000 $quantScript -p -B -t exon -g exon_id -o $quantExonOut/$f.featureCounts.txt $genomeOut/${f}Aligned.sortedByCoord.out.bam  > $quantExonOut/$f.out 2>&1
	else
        	echo "Enter a valid library type. Either 'single' or 'paired'"
        	exit 3
	fi
fi
