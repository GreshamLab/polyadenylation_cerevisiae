#Parameters
GENOME_PATH = "/home/sz4633/polyadenylation_cerevisiae/data/ncbi_dataset/GCA_0001460452/"
GENOME_FASTA = "GCA_0001460452_R64_genomic.fna"
GENOME_GTF = "genomic.gtf"

SAMPLE_PATH = "/scratch/cgsb/gresham/Chris/RAPA_SINGLE_CELL_FASTQ"
SAMPLE = "RAPA1"

WHITELIST = "/scratch/cgsb/gresham/Chris/3M-february-2018.txt"
THREADS = 12

OUTPUT_PATH = "/home/sz4633/polyadenylation_cerevisiae/results/"


#Workflow
rule all:
    input:
        directory(os.path.join(f"{OUTPUT_PATH}", f"{SAMPLE}"))


rule build_genome_index:
    input:
        os.path.join(f"{GENOME_PATH}", f"{GENOME_FASTA}"),
        os.path.join(f"{GENOME_PATH}", f"{GENOME_GTF}")

    output:
        directory(os.path.join(f"{OUTPUT_PATH}", f"star_index"))

    shell:
        """
            mkdir -p {output}

            star_executable/STAR \
             --runThreadN {THREADS} \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input[0]} \
             --sjdbGTFfile {input[1]} \
             --sjdbOverhang 100 \
             --genomeSAindexNbases 10 \
             --outFileNamePrefix genome_index

        """

rule map_fastq_to_genome:
    input:
        os.path.join(f"{OUTPUT_PATH}", f"star_index"),
        os.path.join(f"{SAMPLE_PATH}", f"{SAMPLE}-LIB_R1.fastq.gz"),
        os.path.join(f"{SAMPLE_PATH}", f"{SAMPLE}-LIB_R2.fastq.gz")

    output:
        directory(os.path.join(f"{OUTPUT_PATH}", f"{SAMPLE}/"))

    shell:
        """
            mkdir -p {output}

            star_executable/STAR \
             --runThreadN {THREADS} \
             --genomeDir {input[0]} \
             --readFilesCommand gunzip -c \
             --readFilesIn {input[2]} {input[1]} \
             --soloType CB_UMI_Simple \
             --soloCBwhitelist {WHITELIST} \
             --soloCBstart 1 \
             --soloCBlen 16 \
             --soloUMIstart 17 \
             --soloUMIlen 12 \
             --soloBarcodeReadLength 1 \
             --soloBarcodeMate 0 \
             --soloUMIdedup 1MM_CR \
             --soloUMIfiltering MultiGeneUMI_CR \
             --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
             --soloCellFilter EmptyDrops_CR \
             --outSAMattributes All \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode TranscriptomeSAM GeneCounts \
             --outReadsUnmapped Fastx \
             -- readMapNumber 10000 \
             --outFileNamePrefix {output}/

        """