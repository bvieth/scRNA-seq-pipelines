###################################################################
###########		human		###########################
###################################################################

# download annotation (gtf) and convert to transcript fasta files

# vega
#ftp://ftp.ensembl.org/pub/vega/human/
# reduce to exons
more Homo_sapiens.GRCh38.68_vega.gtf | awk '$3 =="exon"' > Homo_sapiens.GRCh38.68_vega_exon.gtf
# get sequences matching annotation in gtf file
gffread Homo_sapiens.GRCh38.68_vega_exon.gtf -g Homo_sapiens.VEGA68.dna.toplevel.fa -w Homo_sapiens.GRCh38.68_vega_exon.fa

# gencode
#Release 30 (GRCh38.p12): ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz; ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh38.primary_assembly.genome.fa.gz

# reduce to exons
more gencode.v30.primary_assembly.annotation.gtf | awk '$3 =="exon"' > gencode.v30.primary_assembly.annotation_exon.gtf
gffread gencode.v30.primary_assembly.annotation_exon.gtf -g GRCh38.primary_assembly.genome.fa -w gencode.v30.primary_assembly.annotation_exon.fa

# refseq
more Refseqcurated_hg38.gtf | awk '$3 =="exon"' > Refseqcurated_hg38_exon.gtf
gffread Refseqcurated_hg38_exon.gtf -g hg38.fa -w Refseqcurated_hg38_exon.fa
# remove double entries for kallisto and bwa transcriptome mapping
awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' /data/share/htp/powsimR/data/input/annotation/hg38/refseq/Refseqcurated_hg38_exon.fa | awk '!seen[$1]++' | awk -v OFS="\n" '{print $1,$2}' > /data/share/htp/powsimR/data/input/annotation/hg38/refseq/Refseqcurated_hg38_uniqexon.fa

## NOTE gencode and ensembl are the same, continue with gencode

# combine transcript fasta with ERCC fasta
cat Homo_sapiens.GRCh38.68_vega_exon.fa ../../ercc/ercc.fa > hg38_transcriptome_spike_vega.fa
cat gencode.v30.primary_assembly.annotation_exon.fa ../../ercc/ercc.fa > hg38_transcriptome_spike_gencode.fa
cat Refseqcurated_hg38_uniqexon.fa ../../ercc/ercc.fa > hg38_transcriptome_spike_refseq.fa

## transcript / gene id conversion table
awk '{print $10 " " $12}' Homo_sapiens.GRCh38.68_vega_exon.gtf  | sed 's/[";]//g' | sort | uniq > transcript_gene_conv_vega.txt
awk '{print $10 " " $12}' gencode.v30.primary_assembly.annotation_exon.gtf | sed 's/[";]//g' | sort | uniq > transcript_gene_conv_gencode.txt
awk '{print $10 " " $12}' Refseqcurated_hg38_exon.gtf | sed 's/[";]//g' | sort | uniq > transcript_gene_conv_refseq.txt

## dummy ercc transcript gene conversion
awk '{print $1 " " $1}' ercc.gtf > transcript_gene_conv_ercc.txt

## combine transcript conversion with ERCC dummy
cat transcript_gene_conv_vega.txt ../../ercc/transcript_gene_conv_ercc.txt > transcript_gene_conv_spike_vega.txt
cat transcript_gene_conv_gencode.txt ../../ercc/transcript_gene_conv_ercc.txt > transcript_gene_conv_spike_gencode.txt
cat transcript_gene_conv_refseq.txt ../../ercc/transcript_gene_conv_ercc.txt > transcript_gene_conv_spike_refseq.txt

## combine genomes fasta with spike-ins fasta
cat Homo_sapiens.VEGA68.dna.toplevel.fa ../../ercc/ercc.fa > hg38_genome_spike_vega.fa
cat GRCh38.primary_assembly.genome.fa  ../../ercc/ercc.fa > hg38_genome_spike_gencode.fa
cat hg38.fa ../../ercc/ercc.fa > hg38_genome_spike_refseq.fa

## combine genomes gtf with spike-ins gtf
cat Homo_sapiens.GRCh38.68_vega.gtf ../../ercc/ercc.gtf > hg38_genome_spike_vega.gtf
cat gencode.v30.primary_assembly.annotation.gtf ../../ercc/ercc.gtf > hg38_genome_spike_gencode.gtf
cat Refseqcurated_hg38.gtf ../../ercc/ercc.gtf > hg38_genome_spike_refseq.gtf


###################################################################
###########		mouse		###########################
###################################################################

# download annotation (gtf) and convert to transcript fasta files

# vega
#ftp://ftp.sanger.ac.uk/pub/vega/mouse
# reduce to exons
more Mus_musculus.GRCm38.68_vega.gtf | awk '$3 =="exon"' > Mus_musculus.GRCm38.68_vega_exon.gtf
# get sequences matching annotation in gtf file
gffread Mus_musculus.GRCm38.68_vega_exon.gtf -g Mus_musculus.VEGA68.dna.toplevel.fa -w Mus_musculus.GRCm38.68_vega_exon.fa

# gencode
#Release M15 (GRCm38.p5): ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M15/gencode.vM15.primary_assembly.annotation.gtf.gz; ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M15/GRCm38.primary_assembly.genome.fa.gz
# reduce to exons
more gencode.vM15.primary_assembly.annotation.gtf | awk '$3 =="exon"' > gencode.vM15.primary_assembly.annotation_exon.gtf
gffread gencode.vM15.primary_assembly.annotation_exon.gtf -g GRCm38.primary_assembly.genome.fa -w gencode.vM15.primary_assembly_exon.fa

# ensembl
# remove patches
faSomeRecords Mus_musculus.GRCm38.dna.toplevel.fa chromosome.txt Mus_musculus.GRCm38.ensembl.fa
# reduce to exons
more Mus_musculus.GRCm38.90.gtf | awk '$3 =="exon"' > Mus_musculus.GRCm38.90_exon.gtf
gffread Mus_musculus.GRCm38.90_exon.gtf -g Mus_musculus.GRCm38.ensembl.fa -w Mus_musculus.GRCm38.ensembl_exon.fa

# refseq
more Refseqcurated_mm10.gtf | awk '$3 =="exon"' > Refseqcurated_mm10_exon.gtf
gffread Refseqcurated_mm10_exon.gtf -g mm10.fa -w Refseqcurated_mm10_exon.fa

## NOTE gencode and ensembl are the same, continue with gencode


# combine transcript fasta with ERCC fasta
cat Mus_musculus.GRCm38.68_vega_exon.fa ../ercc/ercc.fa > mm10_transcriptome_spike_vega.fa
cat gencode.vM15.primary_assembly_exon.fa ../ercc/ercc.fa > mm10_transcriptome_spike_gencode.fa
cat Mus_musculus.GRCm38.ensembl_exon.fa ../ercc/ercc.fa > mm10_transcriptome_spike_ensembl.fa
cat Refseqcurated_mm10_exon.fa ../ercc/ercc.fa > mm10_transcriptome_spike_refseq.fa
  
## transcript / gene id conversion table
awk '{print $10 " " $12}' Mus_musculus.GRCm38.68_vega_exon.gtf  | sed 's/[";]//g' | sort | uniq > transcript_gene_conv_vega.txt
awk '{print $10 " " $12}' gencode.vM15.primary_assembly.annotation_exon.gtf | sed 's/[";]//g' | sort | uniq > transcript_gene_conv_gencode.txt
awk '{print $10 " " $14}' Mus_musculus.GRCm38.90_exon.gtf | sed 's/[";]//g' | sort | uniq > transcript_gene_conv_ensembl.txt
awk '{print $10 " " $12}' Refseqcurated_mm10_exon.gtf | sed 's/[";]//g' | sort | uniq > transcript_gene_conv_refseq.txt

## dummy ercc transcript gene conversion
awk '{print $1 " " $1}' ercc.gtf > transcript_gene_conv_ercc.txt

## combine transcript conversion with ERCC dummy
cat transcript_gene_conv_ensembl.txt ../ercc/transcript_gene_conv_ercc.txt > transcript_gene_conv_spike_ensembl.txt
cat transcript_gene_conv_vega.txt ../ercc/transcript_gene_conv_ercc.txt > transcript_gene_conv_spike_vega.txt
cat transcript_gene_conv_gencode.txt ../ercc/transcript_gene_conv_ercc.txt > transcript_gene_conv_spike_gencode.txt
cat transcript_gene_conv_refseq.txt ../ercc/transcript_gene_conv_ercc.txt > transcript_gene_conv_spike_refseq.txt

## combine genomes fasta with spike-ins fasta
cat Mus_musculus.VEGA68.dna.toplevel.fa ../ercc/ercc.fa > mm10_genome_spike_vega.fa
cat GRCm38.primary_assembly.genome.fa  ../ercc/ercc.fa > mm10_genome_spike_gencode.fa
cat mm10.fa ../ercc/ercc.fa > mm10_genome_spike_refseq.fa
cat Mus_musculus.GRCm38.ensembl.fa  ../ercc/ercc.fa > mm10_genome_spike_ensembl.fa

## combine genomes gtf with spike-ins gtf
cat Mus_musculus.GRCm38.68_vega.gtf ../ercc/ercc.gtf > mm10_genome_spike_vega.gtf
cat gencode.vM15.primary_assembly.annotation.gtf ../ercc/ercc.gtf > mm10_genome_spike_gencode.gtf
cat Refseqcurated_mm10.gtf ../ercc/ercc.gtf > mm10_genome_spike_refseq.gtf
cat  Mus_musculus.GRCm38.90.gtf ../ercc/ercc.gtf > mm10_genome_spike_ensembl.gtf


