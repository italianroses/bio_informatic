# Manuplate geneexpression data

# load library
library(dplyr)
library(tidyverse)
library(GEOquery)

# read in the data
dat <- read.csv(file = "Desktop/Gene/GSE183947_fpkm.csv")
dim(dat)

# get metadata
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)

# don't know function // Describe
metadata <- pData(phenoData(gse[[1]]))
head(metadata)
view(metadata)

# Filtering Column
metadata.subsets <- select(metadata, c(1,10,11,17))

# save in meatadata no need to create another dataframe
metadata.modified <- metadata %>%
  select(1,10,11,17) %>% # select row
  rename(tissue = characteristics_ch1) %>% # change row name
  rename(meatastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>% # remove tissue: in data
  mutate(meatastasis = gsub("metastasis: ", "", meatastasis))

# reshaping data to long format
View(dat)
dat.long <- dat %>%
  rename(Gene = X) %>%
  gather(key = 'samples', value = 'FPKM', -Gene) # gather(dataframe, key, value, 2:4 // col 2, 3, 4)
  
View(dat.long)
# Join datafram
dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))


# explor data

# Filtering data
head(dat.long)
dat.long %>%
  filter(Gene == 'BRCA1' | Gene == 'BRCA2') %>%
  group_by(Gene, tissue) %>%
  summarize(mead_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) %>%
  arrange(mead_FPKM) # arrange by value // decending order -mead_FPKM


