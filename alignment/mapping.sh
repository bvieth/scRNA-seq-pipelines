#################################################
###########		Human		#########
#################################################	

# 1. indices for all annotations from the transcriptomes with spike-ins
# bwa
bwa index -p /data/share/htp/powsimR/data/input/index/hg38/bwa/gencode/gencode_hg38 /data/share/htp/powsimR/data/input/annotation/hg38/gencode/hg38_transcriptome_spike_gencode.fa
bwa index -p /data/share/htp/powsimR/data/input/index/hg38/bwa/vega/vega_hg38 /data/share/htp/powsimR/data/input/annotation/hg38/vega/hg38_transcriptome_spike_vega.fa
bwa index -p /data/share/htp/powsimR/data/input/index/hg38/bwa/refseq/refseq_hg38 /data/share/htp/powsimR/data/input/annotation/hg38/refseq/hg38_transcriptome_spike_refseq.fa

# kallisto
kallisto index -i /data/share/htp/powsimR/data/input/index/hg38/kallisto/gencode/gencode_hg38 /data/share/htp/powsimR/data/input/annotation/hg38/gencode/hg38_transcriptome_spike_gencode.fa
kallisto index -i /data/share/htp/powsimR/data/input/index/hg38/kallisto/vega/vega_hg38 /data/share/htp/powsimR/data/input/annotation/hg38/vega/hg38_transcriptome_spike_vega.fa
kallisto index -i /data/share/htp/powsimR/data/input/index/hg38/kallisto/refseq/refseq_hg38 /data/share/htp/powsimR/data/input/annotation/hg38/refseq/hg38_transcriptome_spike_refseq.fa

# 2. star indices without annotations but with spike-ins (for zUMIs version zUMI2.4.1) for hg38
# STAR --runMode genomeGenerate --runThreadN 16 --genomeDir mm10_STAR5idx_noGTF --limitGenomeGenerateRAM 111000000000 --genomeFastaFiles mm10.fa
# vega
sbatch --cpus-per-task=12 --mem=50000 -J vegazumi --wrap="/opt/bin/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /data/share/htp/powsimR/data/input/index/hg38/zumi/vega --genomeFastaFiles /data/share/htp/powsimR/data/input/annotation/hg38/vega/hg38_genome_spike_vega.fa --limitGenomeGenerateRAM 111000000000"
# gencode
sbatch --cpus-per-task=14 --mem=50000 -J gencodezumi --wrap="/opt/bin/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /data/share/htp/powsimR/data/input/index/hg38/zumi/gencode --genomeFastaFiles /data/share/htp/powsimR/data/input/annotation/hg38/gencode/hg38_genome_spike_gencode.fa --limitGenomeGenerateRAM 111000000000"
# refseq
sbatch --cpus-per-task=14 --mem=50000 -J refseqzumi --wrap="/opt/bin/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /data/share/htp/powsimR/data/input/index/hg38/zumi/refseq --genomeFastaFiles /data/share/htp/powsimR/data/input/annotation/hg38/refseq/hg38_genome_spike_refseq.fa --limitGenomeGenerateRAM 111000000000"



# 3. map reads to transcriptomes
# see separate slurm scripts in /data/share/htp/powsimR/scripts/slurm/
# slurm_bwa.sh slurm_kallisto_umi.sh slurm_kallisto_noumi.sh

# 4. map reads to genome with STAR for non-UMI methods / with zUMIs for UMI-methods
# see separate slurm scripts in /data/share/htp/powsimR/scripts/slurm/
# slurm_star.sh slurm_zumis.sh slurm_zumis_new.sh

#################################################
###########		Mouse		#########
#################################################	

# 1. indices for all annotations from the transcriptomes with spike-ins
# bwa
bwa index -p /data/share/htp/powsimR/data/input/index/bwa/gencode/gencode_mm10 /data/share/htp/powsimR/data/input/annotation/gencode/mm10_transcriptome_spike_gencode.fa
bwa index -p /data/share/htp/powsimR/data/input/index/bwa/vega/vega_mm10 /data/share/htp/powsimR/data/input/annotation/vega/mm10_transcriptome_spike_vega.fa
bwa index -p /data/share/htp/powsimR/data/input/index/bwa/ensembl/ensembl_mm10 /data/share/htp/powsimR/data/input/annotation/ensembl/mm10_transcriptome_spike_ensembl.fa
bwa index -p /data/share/htp/powsimR/data/input/index/bwa/refseq/refseq_mm10 /data/share/htp/powsimR/data/input/annotation/refseq/mm10_transcriptome_spike_refseq.fa

# bowtie2
bowtie2-build /data/share/htp/powsimR/data/input/annotation/gencode/mm10_transcriptome_spike_gencode.fa /data/share/htp/powsimR/data/input/index/bowtie2/gencode/gencode_mm10
bowtie2-build /data/share/htp/powsimR/data/input/annotation/vega/mm10_transcriptome_spike_vega.fa /data/share/htp/powsimR/data/input/index/bowtie2/vega/vega_mm10
bowtie2-build /data/share/htp/powsimR/data/input/annotation/ensembl/mm10_transcriptome_spike_ensembl.fa /data/share/htp/powsimR/data/input/index/bowtie2/ensembl/ensembl_mm10
bowtie2-build /data/share/htp/powsimR/data/input/annotation/refseq/mm10_transcriptome_spike_refseq.fa /data/share/htp/powsimR/data/input/index/bowtie2/refseq/refseq_mm10

# kallisto
kallisto index -i /data/share/htp/powsimR/data/input/index/kallisto/gencode/gencode_mm10 /data/share/htp/powsimR/data/input/annotation/gencode/mm10_transcriptome_spike_gencode.fa
kallisto index -i /data/share/htp/powsimR/data/input/index/kallisto/vega/vega_mm10 /data/share/htp/powsimR/data/input/annotation/vega/mm10_transcriptome_spike_vega.fa
kallisto index -i /data/share/htp/powsimR/data/input/index/kallisto/ensembl/ensembl_mm10 /data/share/htp/powsimR/data/input/annotation/ensembl/mm10_transcriptome_spike_ensembl.fa
kallisto index -i /data/share/htp/powsimR/data/input/index/kallisto/refseq/refseq_mm10 /data/share/htp/powsimR/data/input/annotation/refseq/mm10_transcriptome_spike_refseq.fa

# 2a. star indices for all annotations with spike-ins for mm10
# vega
sbatch --cpus-per-task=12 --mem=50000 -J vegastar --wrap="/opt/bin/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /data/share/htp/powsimR/data/input/index/star/vega --genomeFastaFiles /data/share/htp/powsimR/data/input/annotation/vega/mm10_genome_spike_vega.fa --sjdbGTFfile /data/share/htp/powsimR/data/input/annotation/vega/mm10_genome_spike_vega.gtf --sjdbOverhang 44"
# gencode
sbatch --cpus-per-task=14 --mem=50000 -J gencodestar --wrap="/opt/bin/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /data/share/htp/powsimR/data/input/index/star/gencode --genomeFastaFiles /data/share/htp/powsimR/data/input/annotation/gencode/mm10_genome_spike_gencode.fa --sjdbGTFfile /data/share/htp/powsimR/data/input/annotation/gencode/mm10_genome_spike_gencode.gtf --sjdbOverhang 44"
# refseq
sbatch --cpus-per-task=14 --mem=50000 -J refseqstar --wrap="/opt/bin/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /data/share/htp/powsimR/data/input/index/star/refseq --genomeFastaFiles /data/share/htp/powsimR/data/input/annotation/refseq/mm10_genome_spike_refseq.fa --sjdbGTFfile /data/share/htp/powsimR/data/input/annotation/refseq/mm10_genome_spike_refseq.gtf --sjdbOverhang 44"
# ensembl
sbatch --cpus-per-task=14 --mem=50000 -J ensemblstar --wrap="/opt/bin/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /data/share/htp/powsimR/data/input/index/star/ensembl --genomeFastaFiles /data/share/htp/powsimR/data/input/annotation/ensembl/mm10_genome_spike_ensembl.fa --sjdbGTFfile /data/share/htp/powsimR/data/input/annotation/ensembl/mm10_genome_spike_ensembl.gtf --sjdbOverhang 44"


# 3a. star indices without annotations but with spike-ins (for zUMIs version zUMI2.4.1) for mm10
# STAR --runMode genomeGenerate --runThreadN 16 --genomeDir mm10_STAR5idx_noGTF --limitGenomeGenerateRAM 111000000000 --genomeFastaFiles mm10.fa
# vega
sbatch --cpus-per-task=12 --mem=50000 -J vegazumi --wrap="/opt/bin/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /data/share/htp/powsimR/data/input/index/mm10/zumi/vega --genomeFastaFiles /data/share/htp/powsimR/data/input/annotation/vega/mm10/mm10_genome_spike_vega.fa --limitGenomeGenerateRAM 111000000000"
# gencode
sbatch --cpus-per-task=14 --mem=50000 -J gencodezumi --wrap="/opt/bin/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /data/share/htp/powsimR/data/input/index/mm10/zumi/gencode --genomeFastaFiles /data/share/htp/powsimR/data/input/annotation/mm10/gencode/mm10_genome_spike_gencode.fa --limitGenomeGenerateRAM 111000000000"
# refseq
sbatch --cpus-per-task=14 --mem=50000 -J refseqzumi --wrap="/opt/bin/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /data/share/htp/powsimR/data/input/index/mm10/zumi/refseq --genomeFastaFiles /data/share/htp/powsimR/data/input/annotation/mm10/refseq/mm10_genome_spike_refseq.fa --limitGenomeGenerateRAM 111000000000"

# 4. map reads to transcriptomes
# see separate slurm scripts in /data/share/htp/powsimR/scripts/slurm/
# slurm_bwa.sh slurm_bowtie2.sh slurm_kallisto_umi.sh slurm_kallisto_noumi.sh

# 5. map reads to genome with STAR for non-UMI methods / with zUMIs for UMI-methods
# see separate slurm scripts in /data/share/htp/powsimR/scripts/slurm/
# slurm_star.sh slurm_zumis.sh slurm_zumis_new.sh




