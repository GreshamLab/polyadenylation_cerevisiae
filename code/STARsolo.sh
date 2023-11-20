#Parameters
GENOME_PATH = "/home/sz4633/polyadenylation_cerevisiae/data/ncbi_dataset/GCA_0001460452/"
GENOME_NAME = "GCA_0001460452_R64_genomic.fna"
ANNOTATION = "genomic.gtf"

FASTQ_PATH = "/scratch/cgsb/gresham/Chris/RAPA_SINGLE_CELL_FASTQ"
OUTPUT_PATH = "/home/sz4633/polyadenylation_cerevisiae/results"

WHITELIST = "/scratch/cgsb/gresham/Chris/3M-february-2018.txt"

THREADS = 12

#Workflow

rule build_genome_index:
    input:
        os.path.join(f"{GENOME_PATH}", f"{GENOME_NAME}"),
        os.path.join(f"{GENOME_PATH}", f"{ANNOTATION}")

    output:
        os.path.join(f"{OUTPUT_PATH}", f"star_index")

    shell:
        """
            STAR --runThreadN {THREADS} \
             --runMode genomeGenerate \
             --genomeDir {GENOME_PATH} \
             --genomeFastaFiles {input[0]} \
             --sjdbGTFfile {input[1]} \
             --sjdbOverhang 100 \
             --outFileNamePrefix genome_index

        """

rule map_fastq_to_genome:
    input:
        os.path.join(f"{OUTPUT_PATH}", f"star_index"),
        os.path.join(f"{FASTQ_PATH}", f"RAPA1-LIB_R1.fastq.gz"),
        os.path.join(f"{FASTQ_PATH}", f"RAPA1-LIB_R2.fastq.gz")

    output:
        os.path.join(f"{OUTPUT_PATH}", f"mapping", f"RAPA1")

    shell:
        """
            STAR --runThreadN {THREADS} \
             --genomeDir {input[0]}\ 
             --readFilesCommand zcat
             --readFilesIn {input[1]} {input[2]} \
             --soloType CB_UMI_Simple \
             --soloCBwhitelist {WHITELIST} \
             --soloCBstart 1 \
             --soloCBlen 16 \
             --soloUMIstart 17 \
             --soloUMIlen\ 12 \
             --soloBarcodeReadLength 1 \
             -- soloBarcodeMate 0 \
             --soloUMIdedup 1MM_CR \
             --soloUMIfiltering MultiGeneUMI_CR \
             --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
             --soloCellFilter EmptyDrops_CR
             --outSAMattributes All \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode TranscriptomeSAM GeneCounts \
             --outReadsUnmapped Fastx \
             --outFileNamePrefix {OUTPUT_PATH}

        """