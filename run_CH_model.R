#######################################################################
# de novo enrichment analysis by CH model
# CH model was firstly described in O'Roak et al. 2012 (PMID:23160955),
# and further developed and applied in Wang et al. 2016 (PMID:27824329) 
# and Coe et al. 2019 (PMID:30559488). Refer above for more details.
#######################################################################

# set up working directory
setwd("/path/to/directory")

# load packages; install first if not
# install.packages("tidyverse")
library(tidyverse)

# load denovop function # FDR correction default as FALSE
source("denovopcallable_default_FALSE.R")

# The CH model gene-level null p-values for LGD and all MISsense
priors <- read.table("CH_Model_null_pvalue_LGD_MIS.txt")

# The CH model gene-level null p-values for LGD and MIS30 only
caddpriors <- read.table("CH_Model_null_pvalue_LGD_MIS30.txt")

# sample size (Below just an example number)
sample_size <- 17426 #  give the total number of probands (or siblings if test for siblings) in analysis

# run CH-model for LGD and MIS
LGD_MIS <- read.table("denovo_LGD_MIS.txt", row.names= 1)
variants <- round(sample_size * 1.8, 0)
CH_model_LGD_MIS <- denovop(genes=LGD_MIS, variants=variants, priors=priors)

# run Ch-model for LGD and MIS30
LGD_MIS30 <- read.table("denovo_LGD_MIS30.txt", row.names= 1)
variants <- round(sample_size * 1.8, 0)
CH_model_LGD_MIS30 <- denovop(genes=LGD_MIS30, variants=variants, priors=caddpriors)

# MERGE P VALUES: LGD, MIS, MIS30
CH_model_LGD_MIS$gene <- row.names(CH_model_LGD_MIS)
CH_model_LGD_MIS30$gene <- row.names(CH_model_LGD_MIS30)

CH_model_combine <- merge(CH_model_LGD_MIS, CH_model_LGD_MIS30,by="gene", all=TRUE)
colnames(CH_model_combine) <- c("gene","dnLGD","dnMIS","dnLGD_pValue","dnMIS_pValue","dnALT_pValue",
                                "dnLGD.2","dnMIS30","dnLGD_pValue.2","dnMIS30_pValue","dnLGD_dnMIS30_pValue")
CH_model_combine_output <- CH_model_combine[,c(1:5,8,10)]

## Multiple test correction of p-Values for all 18946 genes in the CH model using method Benjamini & Hochberg (1995) ("BH" or its alias "fdr")
## NOTE!!! The q-values should be corrected for all genes with DNM in analysis,
## This example code only provided DNM counts for 125 genes as an example.
CH_model_combine_output$dnLGD_qValue <- p.adjust(CH_model_combine_output$dnLGD_pValue,method = "BH", n = 18946)
CH_model_combine_output$dnMIS_qValue <- p.adjust(CH_model_combine_output$dnMIS_pValue,method = "BH", n = 18946)
CH_model_combine_output$dnMIS30_qValue <- p.adjust(CH_model_combine_output$dnMIS30_pValue,method = "BH", n = 18946)

# write to output
CH_model_combine_output <- CH_model_combine_output %>% select(gene, dnLGD, dnMIS, dnMIS30, 
                                                              dnLGD_pValue, dnLGD_qValue,
                                                              dnMIS_pValue, dnMIS_qValue, 
                                                              dnMIS30_pValue, dnMIS30_qValue)

write.csv(CH_model_combine_output,"DNM_CH_model_p_q.csv",row.names = F)

