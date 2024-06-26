#Load packages
library(tidyverse)

#make dir where to save files
if (!dir.exists(snakemake@params[[1]])) {
  dir.create(snakemake@params[[1]])
}

#import data
rawdata <- read.delim(snakemake@input[[1]],
                      header = F,
                      col.names = c("chromosome", 
                                    "peak_start",
                                    "peak_end",
                                    "peak_name",
                                    "gene_name",
                                    "position",
                                    "depth"
                                    )
                      )

df <- rawdata %>%
  mutate(peak_name = as.factor(peak_name)) %>%
  group_by(gene_name) %>%
  filter(depth >= (max(depth)*0.25)) %>% #remove low coverage
  ungroup() %>%
  group_by(chromosome, peak_start, peak_end, peak_name, gene_name) %>%
  summarise(peak_area = sum(depth)) %>% #calculate peak area
  ungroup() %>%
  group_by(gene_name) %>% 
  filter(n() > 1) #remove genes with 1 peak only
  
print(paste("Genes with multiple peaks BEFORE removing low coverage:", 
            length(unique(rawdata$gene_name))
            )
      )

print(paste("Genes with multiple peaks AFTER filtering:", 
            length(unique(df$gene_name))
            )
      )

write_delim(df,
            snakemake@output[[1]], 
            col_names = T,
            delim="\t") 

