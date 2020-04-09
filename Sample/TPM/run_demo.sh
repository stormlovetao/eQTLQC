wget -P ./fastq ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188040/ERR188040_1.fastq.gz
wget -P ./fastq ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188040/ERR188040_2.fastq.gz
wget -P ./fastq ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188231/ERR188231_1.fastq.gz
wget -P ./fastq ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188231/ERR188231_2.fastq.gz
wget -P ./fastq ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188428/ERR188428_1.fastq.gz
wget -P ./fastq ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188428/ERR188428_2.fastq.gz

wget -P ./ref ftp://ftp.ensembl.org/pub/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -P ./ref ftp://ftp.ensembl.org/pub/release-83/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz

gunzip ./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > ./ref/
gunzip ./ref/Homo_sapiens.GRCh38.97.gtf.gz > ./ref/

cd ./fastq/
for gz in *gz; do gunzip $gz; done
cd ../

#create reference
../src/rsem-1.2.25/rsem-prepare-reference --gtf ./ref/Homo_sapiens.GRCh38.97.gtf --bowtie2 --bowtie2-path ../src/bowtie2-2.2.6/ Homo_sapiens.GRCh38.dna.primary_assembly.fa ./ref/human_ensembl

python3 Tool.py







