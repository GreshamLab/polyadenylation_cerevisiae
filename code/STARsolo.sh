#Parameters
GENOME_PATH = "/home/sz4633/polyadenylation_cerevisiae/data/ncbi_dataset/GCA_0001460452/"
GENOME_FASTA = "GCA_0001460452_R64_genomic.fna"
GENOME_GTF = "genomic.gtf"

FASTQ_PATH = "/scratch/cgsb/gresham/Chris/RAPA_SINGLE_CELL_FASTQ"
FASTQ_FILE = ["RAPA1", "RAPA2", "RAPA3", "RAPA4", 
              "RAPA5", "RAPA6", "RAPA7", "RAPA8",
              "RAPA_REP2_1", "RAPA_REP2_2", "RAPA_REP2_3", "RAPA_REP2_4",
              "RAPA_REP2_5", "RAPA_REP2_6", "RAPA_REP2_7", "RAPA_REP2_8"
            ]

WHITELIST = "/scratch/cgsb/gresham/Chris/3M-february-2018.txt"
THREADS = 16

OUTPUT_PATH = "/scratch/sz4633/polyadenylation_cerevisiae/results/"
TMP_DIR = "/scratch/sz4633/polyadenylation_cerevisiae/tmp/"
STAR_PATH = "/home/sz4633/polyadenylation_cerevisiae/code/star_executable"
CODE_FOLDER = "/home/sz4633/polyadenylation_cerevisiae/code"

#Workflow
rule all:
    input:
        expand(os.path.join(f"{OUTPUT_PATH}", "{sample}"), sample = FASTQ_FILE),
        expand(os.path.join(f"{OUTPUT_PATH}", "{sample}/Aligned.out.bam"), sample = FASTQ_FILE),
        expand(os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted.bam"), sample = FASTQ_FILE),
        expand(os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted.bam.bai"), sample = FASTQ_FILE),
        expand(os.path.join(f"{OUTPUT_PATH}", "macs3/{sample}_peaks.xls"), sample = FASTQ_FILE)

    threads: THREADS

rule build_genome_index:
    input:
        os.path.join(f"{GENOME_PATH}", f"{GENOME_FASTA}"),
        os.path.join(f"{GENOME_PATH}", f"{GENOME_GTF}")

    output:
        directory(os.path.join(f"{OUTPUT_PATH}", f"star_index"))

    threads: THREADS

    shell:
        """
            mkdir -p {output}

            {STAR_PATH}/STAR \
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
        os.path.join(f"{FASTQ_PATH}","{sample}-LIB_R1.fastq.gz"),
        os.path.join(f"{FASTQ_PATH}","{sample}-LIB_R2.fastq.gz")

    output:
        directory(os.path.join(f"{OUTPUT_PATH}", "{sample}", "")),
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Aligned.out.bam")
    
    threads: THREADS
    
    shell:
        """
            mkdir -p {output[0]} &&
            cd {output[0]} &&

            {STAR_PATH}/STAR \
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
             --outSAMtype BAM Unsorted \
             --quantMode TranscriptomeSAM GeneCounts \
             --outReadsUnmapped Fastx &&

            cd {CODE_FOLDER}

        """

rule sort_bam_files:
    input:
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Aligned.out.bam")

    output:
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted.bam"),
        temporary(directory(os.path.join(f"{TMP_DIR}", "{sample}")))

    threads: THREADS

    shell:
        """
        mkdir -p {output[1]} &&
        
        samtools sort {input} -o {output[0]} -T {output[1]} -@ {THREADS}

        """

rule index_bam_files:
    input:
        os.path.join(f"{OUTPUT_PATH}", "{sample}", "")

    output:
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted.bam.bai")

    threads: THREADS

    shell:
        """
            cd {input[0]} &&

            samtools index Sorted.bam -@{THREADS} &&

            cd {CODE_FOLDER}

        """

rule find_peaks:
    input:
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted.bam")

    output:
        os.path.join(f"{OUTPUT_PATH}", "macs3/{sample}_peaks.xls")
    
    params:
        outdir = os.path.join(f"{OUTPUT_PATH}", "macs3", "")

    threads: THREADS

    shell:
        """
            macs3 callpeak -t {input} \
                --name {wildcards.sample} \
                --gsize 1.2e7 \
                --nomodel \
                --shift -75 \
                --extsize 150 \
                --outdir {params.outdir}

        """
