#Parameters
GENOME_PATH = "/home/sz4633/polyadenylation_cerevisiae/data/genome_gtf/"
GENOME_FASTA = "Saccharomyces_cerevisiae.R64-1-1.Marker.dna.toplevel.fa"
GENOME_GTF = "Saccharomyces_cerevisiae.R64-1-1.CLEAN.gtf"

FASTQ_PATH = "/scratch/cgsb/gresham/Chris/RAPA_SINGLE_CELL_FASTQ"
FASTQ_FILE = ["RAPA1", "RAPA2", "RAPA3", "RAPA4", 
              "RAPA5", "RAPA6", "RAPA7", "RAPA8",
              "RAPA_REP2_1", "RAPA_REP2_2", "RAPA_REP2_3", "RAPA_REP2_4",
              "RAPA_REP2_5", "RAPA_REP2_6", "RAPA_REP2_7", "RAPA_REP2_8"
            ]

WHITELIST = "/scratch/cgsb/gresham/Chris/3M-february-2018.txt"
INTERSECT_FILE = "/home/sz4633/polyadenylation_cerevisiae/data/Saccharomyces_cerevisiae.R64-1-1.CLEAN.bed"
THREADS = 16

OUTPUT_PATH = "/scratch/sz4633/polyadenylation_cerevisiae/"
TMP_DIR = "/scratch/sz4633/polyadenylation_cerevisiae/tmp/"
STAR_PATH = "/home/sz4633/polyadenylation_cerevisiae/code/star_executable"
CODE_FOLDER = "/home/sz4633/polyadenylation_cerevisiae/code"

#Workflow
rule all:
    input:
        expand(os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted.bam.bai"), sample = FASTQ_FILE),
        os.path.join(f"{CODE_FOLDER}", "check_peaks.html"),
        expand(os.path.join(f"{OUTPUT_PATH}", "quantify_peaks/{sample}"), sample = FASTQ_FILE)


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
             --genomeSAindexNbases 10 `#calculated for small genomes as log2(GenomeLength)/2 - 1`\
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
             --soloCBstart 1 `#barcode start`\
             --soloCBlen 16 `#barcode length`\
             --soloUMIstart 17 `#UMI start`\
             --soloUMIlen 12 `#UMI length`\
             --soloBarcodeReadLength 1 `#length of the barcode read, 1: equal to sum of soloCBlen+soloUMIlen`\
             --soloBarcodeMate 0 `#which read mate contains the barcode: 0 for barcode sequence is on separate read, which should always be the last file in the --readFilesIn listed`\
             --soloUMIdedup 1MM_CR `#type of UMI deduplication (collapsing) algorithm: 1MM_CR for CellRanger2-4 algorithm for 1MM UMI collapsing`\
             --soloUMIfiltering MultiGeneUMI_CR `#remove UMIs with N and homopolymers and remove lower-count UMIs that map to more than one gene`\
             --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts `#option for matching the Cell Barcodes to the WhiteList`\
             --soloCellFilter EmptyDrops_CR `#advanced filtering based on the EmptyDrop algorithm developed by DOI 10.1186/s13059-019-1662-y`\
             --outSAMattributes All \
             --outSAMtype BAM Unsorted `#sorting with star is inefficient, so I sort afterwards with samtools`\
             --quantMode TranscriptomeSAM GeneCounts \
             --outReadsUnmapped Fastx &&

            cd {CODE_FOLDER}

        """

rule sort_bam_files:
#samtools works faster than star for sorting
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

rule index_bam_files_for_IGV:
#index the bam file for viewing in IGV
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
#use MACS3 to find peaks - looks at all the reads present in the sample and finds where peaks are
    input:
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted.bam")

    output:
        os.path.join(f"{OUTPUT_PATH}", "peaks_macs3/{sample}_peaks.narrowPeak")
    
    params:
        outdir = os.path.join(f"{OUTPUT_PATH}", "peaks_macs3", "")

    threads: THREADS

    shell:
        """
            macs3 callpeak -t {input} \
                --name {wildcards.sample} \
                --gsize 1.2e7 `#mappable genome size, defined as the genome size which can be sequenced.` \
                --nomodel `#turn off model construction`\
                --shift -75 `#extend reads from 5’ to 3’ `\
                --extsize 150 `#extend reads in 5’->3’ direction to fix-sized fragments`\
                --max-gap 10 `#max gap between regions`\
                --outdir {params.outdir}

        """

rule bedtools_intersect:
#intersect the peak location (from macs3 bed file) with the gene name in INTERSECT_FILE
    input:
        os.path.join(f"{OUTPUT_PATH}", "peaks_macs3/{sample}_peaks.narrowPeak"),
        os.path.join(f"{INTERSECT_FILE}")

    output:
        os.path.join(f"{OUTPUT_PATH}", "bedtools/{sample}_intersect")

    params:
        outdir = os.path.join(f"{OUTPUT_PATH}", "bedtools", "")
    
    threads: THREADS

    shell:
        """
            cd {params.outdir}

            bedtools intersect -a {input[0]} -b {input[1]} > {output} -wa -wb

            cd {CODE_FOLDER}
        """

rule check_peaks:
#this rule runs the rmd script to filter out and refine which peaks to retain for further analysis
#it produces an Rmarkdown document in which each filtering step is described and justified
    input:
        os.path.join(f"{OUTPUT_PATH}", "bedtools/")

    output:
        os.path.join(f"{CODE_FOLDER}", "check_peaks.html")

    threads: THREADS

    script:
        "check_peaks.Rmd"


rule count_reads_in_bed_file:
#from the bed file in which I have filtered peaks, it counts the reads in the bam alignemnt file
    input:
        os.path.join(f"{OUTPUT_PATH}", "filtered_peaks/{sample}.bed"),
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted.bam")

    output:
        os.path.join(f"{OUTPUT_PATH}", "quantify_peaks/{sample}")

    params:
        outdir = os.path.join(f"{OUTPUT_PATH}", "quantify_peaks", "")

    threads: THREADS

    shell:
        """
            cd {params.outdir}

            bedtools intersect -a {input[0]} -b {input[1]} -c > {wildcards.sample}.bed

            cd {CODE_FOLDER}

        """