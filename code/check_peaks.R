#Load packages
library(tidyverse)

print("Now filtering peaks...")

#make dir where to save files
if (!dir.exists(snakemake@params[[1]])) {
  dir.create(snakemake@params[[1]])
  }
      
#import peak file
rawdata <- read.csv2(snakemake@input[[1]],
                     sep = "\t",
                     header = F,
                     col.names = c("chromosome",
                                   "peak_start",
                                   "peak_end",
                                   "peak_name",
                                   "score",
                                   "strand_macs",
                                   "Fold_Change_at_peak_summit",
                                   "log10pvalue_at_peak_summit",
                                   "log10qvalue_at_peak_summit",
                                   "relative_summit_position_to_peak_start",
                                   "chromosome_position",
                                   "gene_start",
                                   "gene_end",
                                   "gene_name",
                                   "tss",
                                   "strand"
                                   )
                     )

#Load GTF and wrangle 
gtf <- read.table(snakemake@input[[2]],
                  header = FALSE,
                  col.names = c("sequence", 
                                "source", 
                                "feature", 
                                "start", 
                                "end", 
                                "score", 
                                "strand", 
                                "frame", 
                                "gene_name"),
                  sep = '\t') %>% #import data
  mutate(gene_name = str_replace(gene_name, ";.*", ""), #fix gene names
         gene_name = str_replace(gene_name, "gene_id ", "")) %>% #fix gene names
  filter(feature %in% c("transcript", 
                        "start_codon", 
                        "stop_codon")) %>%
  select(-c("sequence", 
            "source", 
            "score", 
            "strand", 
            "frame"))

#Some genes have duplicated start codons that mess up the pivoting. Remove them
duplicated <- gtf %>%
  group_by(gene_name, feature) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1L)

print("These genes have duplicated start codons and will be removed from GTF file")
print(duplicated)

#Pivot the gtf file
gtf <- gtf %>% 
  filter(! gene_name %in% duplicated$gene_name) %>%
  drop_na() %>%
  pivot_wider(names_from = c("feature"),
              values_from = c("start", "end")) %>%
  select(-c("end_start_codon", 
            "end_stop_codon"))

colnames(gtf) <- c("gene_name", 
                   "start_transcript",
                   "start_codon",
                   "stop_codon",
                   "end_transcript")

rm(duplicated)

#Merge peaks_list dataframes with the gtf file. 
peaks <- rawdata %>% 
  select(-c("chromosome_position",
            "strand_macs",
            "gene_start",
            "gene_end",
            "tss")
         ) %>%
  merge(gtf,
        by.x = "gene_name",
        by.y = "gene_name") %>% #remove non important cols
  relocate(gene_name, .after = peak_name)


# Calculate important statistics
  #Transcript length
  #abs_peak_summit: absolute peak position - peak position on chromosome
  #Peak width: info needed to filter out peaks that are too broad
  #Peak to transcript start
  #Peak to transcript end
  #Peak to start codon
  #Peak to end codon

peaks <- peaks %>% 
  mutate(transcript_length = end_transcript - start_transcript,
         abs_peak_summit = peak_start + relative_summit_position_to_peak_start,
         peak_width = peak_end - peak_start, 
         peak_to_transcript_start =
           ifelse(strand == "+", #calculation is different depending on strand
                  abs_peak_summit - start_transcript,
                  end_transcript - abs_peak_summit),
         peak_to_transcript_end = 
           ifelse(strand == "+", 
                  end_transcript - abs_peak_summit,
                  abs_peak_summit - start_transcript),
         peak_to_start_codon =
           ifelse(strand == "+",
                  abs_peak_summit - start_codon,
                  start_codon - abs_peak_summit),
         peak_to_stop_codon = 
           ifelse(strand == "+",
                  stop_codon - abs_peak_summit,
                  abs_peak_summit - stop_codon))

# Check distributions

## Distance to start/stop of transcript/codon
png(snakemake@output[[1]], 
    units = "in", 
    width = 8, 
    height = 4, 
    res = 300)
peaks %>%
      select(c("gene_name",
               "peak_width",
               "peak_to_transcript_start",
               "peak_to_transcript_end",
               "peak_to_start_codon",
               "peak_to_stop_codon")) %>%
      pivot_longer(!gene_name, 
                   names_to = "variable",
                   values_to = "values") %>%
      ggplot(aes(x = values)) + 
      geom_histogram(bins = 500) +
      xlim(-2000,4000) +
      facet_wrap(~variable,
                 scales = "free_y",
                 ncol = 2)
dev.off()

#How many peaks are assigned to more than one gene?
  #Peaks can be broad and can overlap with multiple genes. When this happens, 
  #a single peak (with unique peak name) will be assigned to multiple genes. 
  #For example peak_1 might be assigned to gene A and gene B. 
  #Here we check how many peaks have been assigned to multiple genes. 

tmp <- as.data.frame(table(peaks$peak_name))
tmp <- as.data.frame(table(tmp$Freq))
colnames(tmp) <- c("genes_assigned_to_peak", 
                   "number_of_peaks")
print("How many peaks are assigned to one or more genes?")
tmp

#Check how many genes have more than one peak
  #One single gene might contain multiple peaks. 
  #This is what we are looking for as it is a sign of different 3'UTR.

tmp <- as.data.frame(table(peaks$gene_name))
tmp <- as.data.frame(table(tmp$Freq))
colnames(tmp) <- c("peaks_per_gene", 
                   "frequency")
print("Genes with one or more peaks assigned")
tmp

#Filter out some peaks

  #Remove peaks whose summit is too close to the TSS. Some peaks are broad 
  #enough that one of their tails is very close to the TSS of a gene. 
  #These peaks are for the stop of the other gene, not for the start of a gene. 
  #However bedtools intersect assign that peak to both genes since the peak 
  #technically overlaps - here remove these types of peaks. 

  #Quantify how many of these instances we have:
  print("Peaks whose summit is too close to the TSS. Removing them:")
  nrow(peaks %>% 
       filter(peak_to_transcript_start < 50))
  #Remove them
  peaks <- peaks %>% 
    filter(peak_to_transcript_start > 50)


  #Remove peaks that are too broad
  #Peaks that have a width above the transcript length will 
  #be filtered out since they are too broad

  #Quantify how many of these instances we have:
  tmp1 <- as.data.frame(nrow(peaks %>% 
                               filter(peak_width > 500)
                             )
                        )
  tmp2 <- as.data.frame(nrow(peaks %>% 
                               filter(peak_width > transcript_length)
                             )
                        )

  tmp <- cbind(tmp1,tmp2)
  colnames(tmp) <- c("Peaks_larger_than_500", "Larger_Than_Transcript_Length")
  print("Peaks that are too broad will be removed:")
  tmp
  
  rm(list=ls(pattern="^tmp"))

  #Remove them 
  peaks <- peaks %>% 
    filter(peak_width < transcript_length)

  #Remove genes with only one peak
  #Now that we filtered out most of the unwanted peaks, 
  #let's find the genes with more than one peak! 
  print("Only keeping genes with more than one peak assigned for next steps")
  peaks <- peaks %>% 
    group_by(gene_name) %>% 
    filter(n()>1)
  

#Save output
  #export bed with minumum amount of cols
  write_delim(peaks[1:5],
              snakemake@output[[2]], 
              col_names = F,
              delim="\t") 
  
  #export full df
  write_delim(peaks,
              snakemake@output[[3]], 
              col_names = T,
              delim="\t") 
