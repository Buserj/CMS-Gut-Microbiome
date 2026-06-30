## Construct Venn Diagram of ASVs
## counts shared ASV across environmental compartments
## plots counts using ggVennDiagram

library(ggVennDiagram)
library(tidyverse)

setwd("/wk/dir")
# specify for dataset selection
experiment <- "all" 

#############
# load data #
#############
# select ASV counts dataset
# ASVs as columns, samples as rows
otu <- read.csv(paste0("./CleanData/otu_",experiment,".csv"),row.names = 1)

# import relevant metadata
metadata <- read.table("Metadata.tsv", header = TRUE)
metadata <- metadata[metadata$sample.id %in% rownames(otu), ]

# select for desired groupings
# env (environmental compartment) and guts (snail microbiome)
metadata <- metadata %>% filter(experiment == "lake_guts" | experiment == "ms1_env") 
otu <- otu[rownames(otu) %in% metadata$sample.id, ]
# remove any ASVs with only 0 after filtering
otu <- otu[, colSums(otu !=0)>0]

# separate out lake of interest
trout_metadata <- metadata %>%
  filter(grepl("trout", treatment, ignore.case = TRUE))
trout_otu <- otu[rownames(otu) %in% trout_metadata$sample.id, ]
trout_otu <- trout_otu[, colSums(trout_otu !=0)>0]

# diagram
# trout
trout_asvgroups <- trout_otu %>%
  as.data.frame() %>%
  mutate(treatment = trout_metadata$treatment) %>%
  group_by(treatment) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("treatment")
trout_asvcounts <- apply(trout_asvgroups > 0, 1, function(x) names(which(x)))

names(trout_asvcounts) <- c("Periphyton","Sediment" ,"Water Column","Snail Gut")
ggVennDiagram(trout_asvcounts, label_alpha = 0) +
  scale_fill_gradient(low="steelblue", high = "white", name = "ASV Count")+
  ggtitle("Trout Lake")
