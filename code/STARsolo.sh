#Parameters
GENOME_PATH = "/home/sz4633/polyadenylation_cerevisiae/data/genome_gtf/"
GENOME_FASTA = "Saccharomyces_cerevisiae.R64-1-1.Marker.dna.toplevel.fa"
GENOME_SORTED = "Saccharomyces_cerevisiae.R64-1-1.Marker.dna.txt"
GENOME_GTF = "Saccharomyces_cerevisiae.R64-1-1.CLEAN.gtf"

FASTQ_PATH = "/scratch/cgsb/gresham/Chris/RAPA_SINGLE_CELL_FASTQ"
FASTQ_FILE = ["RAPA1" , "RAPA2"] #, "RAPA3", "RAPA4", 
              #"RAPA5", "RAPA6", "RAPA7", "RAPA8",
              #"RAPA_REP2_1", "RAPA_REP2_2", "RAPA_REP2_3", "RAPA_REP2_4",
              #"RAPA_REP2_5", "RAPA_REP2_6", "RAPA_REP2_7", "RAPA_REP2_8"
            #]
orientation = ["positive" , "negative"]

WHITELIST = "/scratch/cgsb/gresham/Chris/3M-february-2018.txt"
INTERSECT_FILE = "/home/sz4633/polyadenylation_cerevisiae/data/Saccharomyces_cerevisiae.R64-1-1.CLEAN.bed"
THREADS = 16

OUTPUT_PATH = "/scratch/sz4633/polyadenylation_cerevisiae_split_strands/"
TMP_DIR = "/scratch/sz4633/polyadenylation_cerevisiae_split_strands/tmp/"
STAR_PATH = "/home/sz4633/polyadenylation_cerevisiae/code/star_executable"
CODE_FOLDER = "/home/sz4633/polyadenylation_cerevisiae/code"

#Workflow
rule all:
    input:
        expand(os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted_{strand}.bam.bai"), sample = FASTQ_FILE, strand = orientation),
        #expand(os.path.join(f"{OUTPUT_PATH}", "peaks_macs3/{sample}_{strand}_peaks.narrowPeak"), sample = FASTQ_FILE, strand = strand)
        #expand(os.path.join(f"{OUTPUT_PATH}", "peaks_bedtools_intersect/{sample}_{strand}_intersect"), sample = FASTQ_FILE, strand = strand)
        #expand(os.path.join(f"{OUTPUT_PATH}", "peaks_filtered/{sample}_{strand}.tsv"), sample = FASTQ_FILE, strand = strand)
        #expand(os.path.join(f"{OUTPUT_PATH}", "peaks_bedtools_intersect_sorted/{sample}_{strand}_sorted"), sample = FASTQ_FILE, strand = strand)
        #expand(os.path.join(f"{OUTPUT_PATH}", "peaks_seq_depth/{sample}_{strand}"), sample = FASTQ_FILE, strand = strand)
        expand(os.path.join(f"{OUTPUT_PATH}", "peaks_area/{sample}_{strand}_PeakArea.tsv"), sample = FASTQ_FILE, strand = orientation)

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
        genome_index = os.path.join(OUTPUT_PATH, "star_index"),
        fastq_r1 = os.path.join(FASTQ_PATH,"{sample}-LIB_R1.fastq.gz"),
        fastq_r2 = os.path.join(FASTQ_PATH,"{sample}-LIB_R2.fastq.gz")

    output:
        directory(os.path.join(f"{OUTPUT_PATH}", "{sample}")),
        bam_file = os.path.join(OUTPUT_PATH, "{sample}/Aligned.out.bam")
    
    threads: THREADS
    
    shell:
        """
            mkdir -p {output[0]} &&
            cd {output[0]} &&

            {STAR_PATH}/STAR \
             --runThreadN {THREADS} \
             --genomeDir {input.genome_index} \
             --readFilesCommand gunzip -c \
             --readFilesIn {input.fastq_r2} {input.fastq_r1} \
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

rule split_strands: #split reads that align to positive and negative strand
    input:
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Aligned.out.bam")

    output:
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Aligned_{strand}.out.bam"),

    threads: THREADS
    
    shell:
        """
        if [ "{wildcards.strand}" = "negative" ]
        then
            samtools view -f 16 -o {output} {input}
        else
            samtools view -F 16 -o {output} {input}
        fi

        """

rule sort_bam_files: #samtools works faster than STAR for sorting
    input:
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Aligned_{strand}.out.bam")

    output:
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted_{strand}.bam"),
        temporary(directory(os.path.join(f"{TMP_DIR}", "{sample}_{strand}")))

    threads: THREADS

    shell:
        """
        mkdir -p {output[1]} &&
        
        samtools sort {input} -o {output[0]} -T {output[1]} -@ {THREADS}

        """

rule index_bam_files_for_IGV: #index the sorted.bam file for IGV
    input:
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted_{strand}.bam")

    output:
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted_{strand}.bam.bai")
    
    threads: THREADS

    params:
        outdir = os.path.join(f"{OUTPUT_PATH}", "{sample}/")
    
    shell:
        """
            cd {params.outdir} &&

            samtools index {input} -@{THREADS} &&

            cd {CODE_FOLDER}

        """

rule find_peaks: #looks at reads in bam file and finds peaks
    input:
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted_{strand}.bam"),

    output:
        os.path.join(f"{OUTPUT_PATH}", "peaks_macs3/{sample}_{strand}_peaks.narrowPeak")
    
    params:
        outdir = os.path.join(f"{OUTPUT_PATH}", "peaks_macs3", "")

    threads: THREADS

    shell:
        """
            macs3 callpeak -t {input} \
                --name {wildcards.sample}_{wildcards.strand} \
                --gsize 1.2e7 `#mappable genome size, defined as the genome size which can be sequenced.` \
                --nomodel `#turn off model construction`\
                --shift -75 `#extend reads from 5’ to 3’ `\
                --extsize 150 `#extend reads in 5’->3’ direction to fix-sized fragments`\
                --max-gap 10 `#max gap between regions`\
                --outdir {params.outdir}

        """

rule bedtools_intersect: #intersect peak location with the gene name
    input:
        os.path.join(f"{OUTPUT_PATH}", "peaks_macs3/{sample}_{strand}_peaks.narrowPeak"),
        os.path.join(f"{INTERSECT_FILE}")

    output:
        os.path.join(f"{OUTPUT_PATH}", "peaks_bedtools_intersect/{sample}_{strand}_intersect")

    params:
        outdir = os.path.join(f"{OUTPUT_PATH}", "peaks_bedtools_intersect", "")
    
    threads: THREADS

    shell:
        """
            cd {params.outdir}

            bedtools intersect -a {input[0]} -b {input[1]} > {output} -wa -wb

            cd {CODE_FOLDER}
        """

rule filter_interesting_peaks: #filter which peaks to retain for further analysis
    input:
        os.path.join(f"{OUTPUT_PATH}", "peaks_bedtools_intersect/{sample}_{strand}_intersect"),
        os.path.join(f"{GENOME_PATH}", f"{GENOME_GTF}")

    output:
        os.path.join(f"{OUTPUT_PATH}", "peaks_filtered/{sample}_{strand}_plot.png"),
        os.path.join(f"{OUTPUT_PATH}", "peaks_filtered/{sample}_{strand}.tsv"),
        os.path.join(f"{OUTPUT_PATH}", "peaks_filtered/{sample}_{strand}_full.tsv")
    
    params:
        os.path.join(f"{OUTPUT_PATH}", "peaks_filtered/", "")

    threads: THREADS

    script:
        "check_peaks.R"

rule sort_filtered_peaks: #sort bed file in the same way genome and bam files are
    input:
        os.path.join(f"{OUTPUT_PATH}", "peaks_filtered/{sample}_{strand}.tsv"),
        os.path.join(f"{GENOME_PATH}", f"{GENOME_SORTED}")

    output:
        os.path.join(f"{OUTPUT_PATH}", "peaks_bedtools_intersect_sorted/{sample}_{strand}_sorted")

    threads: THREADS

    shell:
        """
            bedtools sort -i {input[0]} -g {input[1]} > {output}

        """

rule calculate_seq_depth: #calculate sequencing depth, and later use it to trim peaks
    input:
        os.path.join(f"{GENOME_PATH}", f"{GENOME_SORTED}"),
        os.path.join(f"{OUTPUT_PATH}", "peaks_bedtools_intersect_sorted/{sample}_{strand}_sorted"),
        os.path.join(f"{OUTPUT_PATH}", "{sample}/Sorted_{strand}.bam")

    output:
        os.path.join(f"{OUTPUT_PATH}", "peaks_seq_depth/{sample}_{strand}")

    threads: THREADS

    shell:
        """

            bedtools coverage -sorted \
                              -d \
                              -split \
                              -g {input[0]} \
                              -a {input[1]} \
                              -b {input[2]} \
                              > {output}

        """

rule trim_peaks: 
    input:
        os.path.join(f"{OUTPUT_PATH}", "peaks_seq_depth/{sample}_{strand}")

    output:
        os.path.join(f"{OUTPUT_PATH}", "peaks_area/{sample}_{strand}_PeakArea.tsv")
    
    params:
        os.path.join(f"{OUTPUT_PATH}", "peaks_area/", "")

    threads: THREADS

    script:
        "trim_peaks.R"

