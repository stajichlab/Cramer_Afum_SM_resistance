#!/usr/bin/env Rscript
library(purrr)
library(tidyverse)

dt <- read_tsv('snpEff/Af293/All/Afum_SM_v1.All.snpEff.matrix.tsv',col_names=TRUE) %>% 
  select(-c(FLANKING, ANN)) %>% mutate_all(function(x) { str_replace(x, "/", "")})
SM12 <- c('SM12_R10', 'SM12_R2', 'SM12_R3', 'SM12_R4', 'SM12_R6', 'SM12_R7')
SM6 <- c('SM6_R1', 'SM6_R2',	'SM6_R3',	'SM6_R4',	'SM6_R6')

m <- dt %>% filter(TYPE != 'intergenic') %>% select(!c(IMPACT,CHANGEPEP,CHANGEDNA,ALT)) %>% 
  pivot_longer(!c(CHROM,POS,TYPE,GENE), names_to="isolate", values_to="genotype") %>%
  separate(isolate, c("Group"),sep="_",extra="drop" ) %>% group_by(CHROM,POS,TYPE,GENE,Group) %>% 
  summarise(genotype = paste(genotype, collapse = ",")) %>% 
  pivot_wider(id_cols=c(CHROM,POS,TYPE,GENE),names_from="Group",values_from="genotype")
write_tsv(m,"variants_grouped.tsv")
