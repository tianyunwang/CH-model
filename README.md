# CH model
> _de novo_ enrichment analysis by CH model

- ### Introduction
Briefly, the CH model estimates the number of expected DNMs by incorporating locus-specific transition, transversion, and indel rates and chimpanzee–human coding sequence divergence and the gene length. The expected mutation rate of 1.8 DNMs per exome was set to the CH model as an upper bound baseline. The CH model was firstly developed and described in O'Roak et al. 2012 (PMID: 23160955), and further developed and applied in Wang et al. 2016 (PMID: 27824329) and Coe et al. 2019 (PMID: 30559488). Also refer above publications for more details.

- ### Genes and null p-values
The CH model contains 18,946 genes in model (hg19), the null p-values are pre-calculated for both de novo likely-gene distruptive (LGD) and missense (MIS) mutations (file `CH_Model_null_pvalue_LGD_MIS.txt`), and also for severe missense mutations based on a combined annotation dependent depletion (CADD) score (Kircher et al., 2014 - PMID: 24487276; Rentzsch et al., 2018 - PMID: 30371827) greater than 30 (v1.3), which ranking these mutations among the top 0.1% for the most severe predicted effect (MIS30) (file `CH_Model_null_pvalue_LGD_MIS30.txt`).

- ### How to run
Prepare your input files of _de novo_ mutation count (LGD + MIS mutation count as as in example file `denovo_LGD_MIS.txt`, and LGD + MIS30 mutation count as in example file `denovo_LGD_MIS30.txt`); then give the total number of probands (or siblings if test for siblings) in analysis as `sample_size` in `run_CH_model.R` and run.
