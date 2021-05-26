import os
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

workingDirectory = ''           # Set the path to your working directory here
pathToPacbioReads = ''          # Set the path to the fastq file containing the PacBio long reads
pathToIlluminaReads1 = ''       # Set the path to the fastq file containing the Illumina short reads pair 1
pathToIlluminaReads2 = ''       # Set the path to the fastq file containing the Illumina short reads pair 2

os.mkdir(f'{workingDirectory}/longqc/')
os.system(f'cd {workingDirectory}/longqc/; git clone https://github.com/yfukasawa/LongQC.git')
os.system(f'cd {workingDirectory}/longqc/LongQC/minimap2_mod && make extra')

##################################################
##### Quality Control of the long reads
##################################################

os.system(f'python longQC.py sampleqc -x pb-sequel -o out_dir {pathToPacbioReads}')

##################################################
##### Check for read contaminations with Kraken2
##################################################

os.mkdir(f'{workingDirectory}/kraken2_db/')


#make customized kraken2 library by adding three diatoms genomes that are relatively close to the sequenced species
os.mkdir(f'{workingDirectory}/additional_diatoms_genomes/')
os.system(f'wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/405/GCF_000149405.2_ASM14940v2/GCF_000149405.2_ASM14940v2_genomic.fna.gz -O {workingDirectory}/additional_diatoms_genomes/thalassiosira_pseudonana_genome.fna.gz')
os.system(f'gunzip {workingDirectory}/additional_diatoms_genomes/thalassiosira_pseudonana_genome.fna.gz')
os.system(f'wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/AG/NL/AGNL01/AGNL01.1.fsa_nt.gz -O {workingDirectory}/additional_diatoms_genomes/thalassiosira_oceanica_contigs.fna.gz')
os.system(f'gunzip {workingDirectory}/additional_diatoms_genomes/thalassiosira_oceanica_contigs.fna.gz')

# The skeletonema costatum genome assembly was obtained from: Ogura, A. et al. Comparative genome and transcriptome analysis of diatom, Skeletonema costatum, reveals evolution of genes for harmful algal bloom. BMC Genomics 19, 765 (2018).
# The sequence headers had to be adjusted to be processed by kraken2
cache = open(f'{workingDirectory}/additional_diatoms_genomes/skeletonema_costatum_contigs.fna').readlines()
with open(f'{workingDirectory}/additional_diatoms_genomes/skeletonema_costatum_contigs.fna', 'w') as outFile:
    for line in cache:
        if line[0] == '>':
            outFile.write(f'{line.strip()}|kraken:taxid|2843\n')
        else:
            outFile.write(line)


os.system(f'nice kraken2-build --download-taxonomy --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-library archaea --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-library bacteria  --threads 48 --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-library plasmid  --threads 48 --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-library viral  --threads 48 --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-library human  --threads 48 --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-library protozoa  --threads 48 --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --download-library UniVec  --threads 48 --db {workingDirectory}/kraken2_db/')


os.system(f'kraken2-build --add-to-library {workingDirectory}/additional_diatoms_genomes/thalassiosira_pseudonana_genome.fna --db {workingDirectory}/kraken2_db/')
os.system(f'kraken2-build --add-to-library {workingDirectory}/additional_diatoms_genomes/thalassiosira_oceanica_contigs.fna --db {workingDirectory}/kraken2_db/')
os.system(f'kraken2-build --add-to-library {workingDirectory}/additional_diatoms_genomes/skeletonema_costatum_contigs.fna --db {workingDirectory}/kraken2_db/')


os.system(f'nice kraken2-build --build --threads 48  --db {workingDirectory}/kraken2_db/')
os.system(f'nice kraken2-build --clean --threads 48  --db {workingDirectory}/kraken2_db/')

#parse sequencing fastq file to fasta file for kraken2 input
with open(f'{workingDirectory}/pacbio_long_reads.fasta', 'w') as outputFile:
    with open(f'{pathToPacbioReads}') as inputFile:
        count = 0
        for line in inputFile:
            count += 1
            if count%4 in (1,2):
                outputFile.write(line.replace('@','>'))

#execute kraken2 on the long read data
os.system(f'nice kraken2 --db {workingDirectory}/kraken2_db/ --threads 48 --use-names --output {workingDirectory}/kraken2_db/kraken2_output_long_reads --report {workingDirectory}/kraken2_db/kraken2_report_long_reads {workingDirectory}/pacbio_long_reads.fasta')

#execute kraken2 on the short read data
os.system(f'nice kraken2 --db {workingDirectory}/kraken2_db/ --threads 48 --use-names --output {workingDirectory}/kraken2_db/kraken2_output_short_reads_L1 --report {workingDirectory}/kraken2_db/kraken2_report_short_reads_L1_1 {pathToIlluminaReads1}')
os.system(f'nice kraken2 --db {workingDirectory}/kraken2_db/ --threads 48 --use-names --output {workingDirectory}/kraken2_db/kraken2_output_short_reads_L2 --report {workingDirectory}/kraken2_db/kraken2_report_short_reads_L1_2 {pathToIlluminaReads2}')


#analyzed the kraken2 report using the pavian web app: https://fbreitwieser.shinyapps.io/pavian/
#if the web app should not work or not be available, it can be installed relatively easy locally: https://github.com/fbreitwieser/pavian


# get the reads that are classified as Skeletonema costatum, Thalassiosira oceanica, Thalassiosira pseudonana and unclassified 
os.mkdir(f'{workingDirectory}/filtered_reads/')

filteredReads = set()

with open(f'{workingDirectory}/kraken2_db/kraken2_output_long_reads') as infile:
    for line in infile:
        if line[0] == 'U':
            filteredReads.add(line.split()[1])
        elif 'Thalassi' in line.split('\t')[2]:
            filteredReads.add(line.split()[1])
        elif 'Skeletonema' in line.split('\t')[2]:
            filteredReads.add(line.split()[1])

writeSeq = False
with open(f'{workingDirectory}/filtered_reads/pacbio_long_reads.filtered.fastq', 'w') as outputFile:
    for line in open(f'{pathToPacbioReads}'):
        if line[0] == '@':
            writeSeq = False
            if line[1:].strip() in filteredReads:
                writeSeq = True
        if writeSeq:
            outputFile.write(line)
            outputFile.flush()


##################################################
##### Assembly using flye
##################################################

os.mkdir(f'{workingDirectory}/assembly/')
os.system(f'nice flye --pacbio-raw {workingDirectory}/filtered_reads/pacbio_long_reads.filtered.fastq  -g 30m -o {workingDirectory}/assembly/ -t 48')


#Use Kraken2 to check the assembled contigs
os.system(f'nice kraken2 --db {workingDirectory}/kraken2_db/ --threads 48 --use-names --output {workingDirectory}/kraken2_db/kraken2_output_assembly --report {workingDirectory}/kraken2_db/kraken2_report_assembly {workingDirectory}/assembly/assembly.fasta')


#Remove the non-diatom contigs from the assembly
filteredContigs = set()
with open(f'{workingDirectory}/kraken2_db/kraken2_output_assembly') as infile:
    for line in infile:
        if 'Thalassi' in line.split('\t')[2]:
            filteredContigs.add(line.split()[1])
        elif 'Skeletonema' in line.split('\t')[2]:
            filteredContigs.add(line.split()[1])


os.makedirs(f'{workingDirectory}/assembly/cleaned/')
with open(f'{workingDirectory}/assembly/assembly_cleaned.fasta', 'w') as outputFile:
    writeSeq = False
    for line in open(f'{workingDirectory}/assembly/assembly.fasta'):
        if line[0] == '>': 
            if line[1:].strip() in filteredContigs:
                outputFile.write(f'{line.strip()}\n')
                writeSeq = True
            else:
                writeSeq = False
        else:
            if writeSeq:
                outputFile.write(f'{line.strip()}\n')
                outputFile.flush()

##################################################
##### Final Assembly polishing with pilon and the Illumina short reads
##################################################

#remove all non-diatom reads from the short read libraries
filteredReadsL1 = set()
with open(f'{workingDirectory}/kraken2_db/kraken2_output_short_reads_L1') as infile:
    for line in infile:
        if 'Thalassi' in line.split('\t')[2]:
            filteredReadsL1.add(line.split()[1])
        elif 'Skeletonema' in line.split('\t')[2]:
            filteredReadsL1.add(line.split()[1])

writeSeq = False
with open(f'{workingDirectory}/filtered_reads/illumina_short_reads_L1.filtered.fastq', 'w') as outputFile:
    for line in open(f'{pathToIlluminaReads1}'):
        if line[0] == '@':
            writeSeq = False
            if line[1:].split()[0] in filteredReadsL1:
                writeSeq = True
        if writeSeq:
            outputFile.write(line)
            outputFile.flush()


filteredReadsL2 = set()
with open(f'{workingDirectory}/kraken2_db/kraken2_output_short_reads_L2') as infile:
    for line in infile:
        if 'Thalassi' in line.split('\t')[2]:
            filteredReadsL2.add(line.split()[1])
        elif 'Skeletonema' in line.split('\t')[2]:
            filteredReadsL2.add(line.split()[1])

writeSeq = False
with open(f'{workingDirectory}/filtered_reads/illumina_short_reads_L2.filtered.fastq', 'w') as outputFile:
    for line in open(f'{pathToIlluminaReads2}'):
        if line[0] == '@':
            writeSeq = False
            if line[1:].split()[0] in filteredReadsL2:
                writeSeq = True
        if writeSeq:
            outputFile.write(line)
            outputFile.flush()


#build hisat2 index for mapping
os.makedirs(f'{workingDirectory}/assembly/polished/')
os.system(f'hisat2-build -p 48 {workingDirectory}/assembly/assembly_cleaned.fasta {workingDirectory}/assembly/polished/assembly_cleaned')

#short read mapping with hisat2
os.system(f'hisat2 --no-spliced-alignment --summary-file {workingDirectory}/assembly/polished/short_read_mapping_L1_summary.txt --new-summary -p 48 -x {workingDirectory}/assembly/polished/assembly_cleaned -U {workingDirectory}/filtered_reads/illumina_short_reads_L1.filtered.fastq -S {workingDirectory}/assembly/polished/short_read_mapping_L1.sam')
os.system(f'hisat2 --no-spliced-alignment --summary-file {workingDirectory}/assembly/polished/short_read_mapping_L2_summary.txt --new-summary -p 48 -x {workingDirectory}/assembly/polished/assembly_cleaned -U {workingDirectory}/filtered_reads/illumina_short_reads_L2.filtered.fastq -S {workingDirectory}/assembly/polished/short_read_mapping_L2.sam')

#sam to bam conversion, merging, sorting and indexing
os.system(f'samtools view -@ 48 -b {workingDirectory}/assembly/polished/short_read_mapping_L1.sam -o {workingDirectory}/assembly/polished/short_read_mapping_L1.bam')
os.system(f'samtools view -@ 48 -b {workingDirectory}/assembly/polished/short_read_mapping_L2.sam -o {workingDirectory}/assembly/polished/short_read_mapping_L2.bam')
os.system(f'samtools merge -@ 48 {workingDirectory}/assembly/polished/short_read_mapping_merged.bam {workingDirectory}/assembly/polished/short_read_mapping_L1.bam {workingDirectory}/assembly/polished/short_read_mapping_L2.bam')
os.system(f'samtools sort -@ 48 {workingDirectory}/assembly/polished/short_read_mapping_merged.bam -o {workingDirectory}/assembly/polished/short_read_mapping_merged_sorted.bam')
os.system(f'samtools index {workingDirectory}/assembly/polished/short_read_mapping_merged_sorted.bam')

#assembly polishing with pilon 
os.makedirs(f'{workingDirectory}/assembly/polished/pilon_output/')
os.system('JAVA_TOOL_OPTIONS="-Xmx400G -Xss2560k"')
os.system(f'java -jar pilon-1.23.jar --threads 48 --genome {workingDirectory}/assembly/assembly_cleaned.fasta --unpaired {workingDirectory}/assembly/polished/short_read_mapping_merged_sorted.bam --outdir {workingDirectory}/assembly/polished/pilon_output/ --output final_assembly_polished.fna')


#final quast report
os.system(f'quast -t 48 --fragmented --glimmer -o {workingDirectory}/assembly/polished/pilon_output/ {workingDirectory}/assembly/polished/pilon_output/final_assembly_polished.fna')