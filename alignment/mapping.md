Alignment
================

# Creating Indices

This is for the annotations created in annotation.Rmd for BWA, kallisto
and
STAR.

## MOUSE

### BWA

``` bash
bwa index -p index/mm10/bwa/gencode/gencode_mm10 annotation/mm10/gencode/mm10_transcriptome_spike_gencode.fa
bwa index -p index/mm10/bwa/vega/vega_mm10 annotation/mm10/vega/mm10_transcriptome_spike_vega.fa
bwa index -p index/mm10/bwa/refseq/refseq_mm10 annotation/mm10/refseq/mm10_transcriptome_spike_refseq.fa
```

### kallisto

``` bash
kallisto index -i index/mm10/kallisto/gencode/gencode_mm10 annotation/mm10/gencode/mm10_transcriptome_spike_gencode.fa
kallisto index -i index/mm10/kallisto/vega/vega_mm10 annotation/mm10/vega/mm10_transcriptome_spike_vega.fa
kallisto index -i index/mm10/kallisto/refseq/refseq_mm10 annotation/mm10/refseq/mm10_transcriptome_spike_refseq.fa
```

### STAR

STAR genome index creation for zUMIs without annotation
information.

``` bash
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir index/mm10/zumi/vega --genomeFastaFiles annotation/vega/mm10/mm10_genome_spike_vega.fa --limitGenomeGenerateRAM 111000000000
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir index/mm10/zumi/gencode --genomeFastaFiles annotation/mm10/gencode/mm10_genome_spike_gencode.fa --limitGenomeGenerateRAM 111000000000
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir index/mm10/zumi/refseq --genomeFastaFiles annotation/mm10/refseq/mm10_genome_spike_refseq.fa --limitGenomeGenerateRAM 111000000000
```

STAR indices with
annotation.

``` bash
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir index/mm10/star/vega --genomeFastaFiles annotation/mm10/vega/mm10_genome_spike_vega.fa --sjdbGTFfile annotation/mm10vega/mm10_genome_spike_vega.gtf --sjdbOverhang 44
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir index/mm10/star/gencode --genomeFastaFiles annotation/mm10/gencode/mm10_genome_spike_gencode.fa --sjdbGTFfile annotation/mm10/gencode/mm10_genome_spike_gencode.gtf --sjdbOverhang 44
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir index/mm10/star/refseq --genomeFastaFiles annotation/mm10/refseq/mm10_genome_spike_refseq.fa --sjdbGTFfile annotation/mm10/refseq/mm10_genome_spike_refseq.gtf --sjdbOverhang 44
```

## Human

### BWA

``` bash
bwa index -p index/hg38/bwa/gencode/gencode_hg38 annotation/hg38/gencode/hg38_transcriptome_spike_gencode.fa
bwa index -p index/hg38/bwa/vega/vega_hg38 annotation/hg38/vega/hg38_transcriptome_spike_vega.fa
bwa index -p index/hg38/bwa/refseq/refseq_hg38 annotation/hg38/refseq/hg38_transcriptome_spike_refseq.fa
```

### kallisto

``` bash
kallisto index -i index/hg38/kallisto/gencode/gencode_hg38 annotation/hg38/gencode/hg38_transcriptome_spike_gencode.fa
kallisto index -i index/hg38/kallisto/vega/vega_hg38 annotation/hg38/vega/hg38_transcriptome_spike_vega.fa
kallisto index -i index/hg38/kallisto/refseq/refseq_hg38 annotation/hg38/refseq/hg38_transcriptome_spike_refseq.fa
```

### STAR

STAR genome index creation for zUMIs without annotation
information.

``` bash
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir index/hg38/zumi/vega --genomeFastaFiles annotation/hg38/vega/hg38_genome_spike_vega.fa --limitGenomeGenerateRAM 111000000000
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir index/hg38/zumi/gencode --genomeFastaFiles annotation/hg38/gencode/hg38_genome_spike_gencode.fa --limitGenomeGenerateRAM 111000000000
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir index/hg38/zumi/refseq --genomeFastaFiles annotation/hg38/refseq/hg38_genome_spike_refseq.fa --limitGenomeGenerateRAM 111000000000
```

# Transcriptome alignment

## BWA

``` bash
# settings
annotation = "gencode"
species = "mm10"
inpath = "path/to/input.fastq/"
fqfile = "cell.cDNA.reads.fastq" # demultiplexed cDNA fastq files, i.e. per cell barcodes
threads = 2
# alignment commands
bwa aln -t $threads index/$species/bwa/$annotation/{$annotation}_{$genome} $inpath/$fqfile > $fqfile.sai
bwa samse index/$species/bwa/$annotation/{$annotation}_{$genome} $fqfile.sai $inpath/$fqfile > $fqfile.sam
samtools view -bS $fqfile.sam > $fqfile.bam
```

# Pseudoalignment

## kallisto

``` bash
# settings
annotation = "gencode"
genome = "mm10"
outpath = "path/to/create/outputfiles/"
threads = 2

# for UMI methods (SCRB-seq, CEL-seq2, Drop-seq, 10X Genomics)
sdfrag = 200 # SCRB-seq
meanfrag = 500 # SCRB-seq
infiles = "inputfiles.txt" # input text file with 1. column = BC; 2. column = UMI.reads.txt (nucleotide sequence of UMI as a text file); 3. column = path and file name of input cDNA reads in fastq format
kallisto pseudo -i index/$species/kallisto/$annotation/{$annotation}_{$genome} -o $outpath -b $fqpath/$protocol$umi --single --umi -t $threads  -l $meanfrag -s $sdfrag

# for Smart-seq2
sdfrag = 200 # Smart-seq2
meanfrag = 500 # Smart-seq2
infiles = "inputfiles.txt" # input text file with 1. column = BC; 2. column = UMI.reads.txt (nucleotide sequence of UMI as a text file); 3. column = path and file name of input cDNA reads in fastq format
infiles = "inputfiles.txt" # input text file with 1. column = BC; 2. column = path and file name of input cDNA reads in fastq format
kallisto pseudo --index index/$species/kallisto/$annotation/{$annotation}_{$genome} --output-dir $outpath --batch $infiles --single --fragment-length $meanfrag --sd $sdfrag --threads $threads
```

# Genome alignment

## zUMIs

Create a yaml config file (see example in alignment folder) and then run
zUMIs:

``` bash
bash /bin/zUMIs/zUMIs-master.sh -y alignmentexample_scrbseq.yaml
```

For Smart-seq2 data simply leave out UMI definition in the sequence file
definition (yaml).

``` bash
```
