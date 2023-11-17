#Parameters
GENOME_PATH = "/data/ncbi_dataset/GCA_000146045.2/"
GENOME_NAME = "GCA_000146045.2_R64_genomic.fna"
ANNOTATION = "genomics.gtf"

FASTQ_PATH = "/scratch/cgsb/gresham/Chris/RAPA_SINGLE_CELL_FASTQ"
OUTPUT_PATH = "/results/STAR"

WHITELIST = "/scratch/cgsb/gresham/Chris/3M-february-2018.txt"

THREADS = 6

#Workflow

rule build_genome_index:
    input:
        os.path.join(f"{GENOME_PATH}", f"{GENOME_NAME}")
        os.path.join(f"{GENOME_PATH}", f"{ANNOTATION}")

    outout:
        os.path.join(f"{OUTPUT_PATH}", f"genome_index")

    shell:
        """
        STAR --runThreadN {THREADS} \
             --runMode genomeGenerate \
             --genomeDir {GENOME_PATH} \
             --genomeFastaFiles {GENOME_PATH}/{GENOME_NAME} \
             --sjdbGTFfile {GENOME_PATH}/{ANNOTATION} \
             --sjdbOverhang ReadLength-1 \
             --outFileNamePrefix genome_index

        """

rule map_fastq_to_genome
    input
        os.path.join(f"{FASTQ_PATH}", f"RAPA1-LIB_R1.fastq.gz")
        os.path.join(f"{FASTQ_PATH}", f"RAPA1-LIB_R2.fastq.gz")
        os.path.join(f"{OUTPUT_PATH}", f"genome_index")
    output
        os.path.join(f"{OUTPUT_PATH}", f"RAPA1")
    shell:
    """
    STAR --runThreadN {THREADS} \
     --genomeDir {input[2]}\ 
     --readFilesCommand zcat
     --readFilesIn {input[0]} {input[1]} \
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