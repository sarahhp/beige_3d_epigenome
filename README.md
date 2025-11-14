Code and analyzed data files to accompany paper Hazell Pickering et al 2025 "Multimodal epigenetic and enhancer network remodeling shape the transcriptional landscape of human beige adipocytes"

in press, Communications Biology 2025

\[https://www.biorxiv.org/content/10.1101/2025.03.28.645896v1]



**Figure 1 RNA-Seq**

Differential Expression (1E)

* RNA-seq/limma\_DE\_donor6\_endpt.html

GSEA comparison with paired BAT and WAT biopsies from Din at al. 2018 Fig S2

* RNA-seq/din2018\_GSEA\_comparsions.html

Deconvolution using snRNA-seq from Miranda et al. Fig S2 2

* RNA-seq/make\_custom\_sigmatrix\_miranda2025.html
* RNA-seq/epic\_plots.html



**Figure 2 ATAC-seq**
(and chromatin states)

* ATAC-seq/ATAC_pipeline.bash

DAR calling (2A)
* ATAC-seq/ATAC_DARs_calling.Rmd

Chromatin states (Supp 5)

* ChIPseq/

Gene expression in DARs (2F)

* ATAC-seq/gene\_expr\_in\_DARs.html



**Figure 3 SHAP integration**

Normalising RNA-seq reads in promoters

* RNA-seq/promoter\_counts\_processing.html

Python notebooks

* SHAP\_integration/py\_nb/

Parameter testing and Random forest comparsion (Supp S7)

* plots\_hyperparameter\_search.html
* random\_forest\_A6H\_7.html

Optimised lightGBM model with SHAP interpretations (3B-D)

* A6H\_7\_features\_opt\_params.html
* waterfall\_figures $$$ better figs elewhere? (3I,J)

Clustering SHAP values (3E-H)

* SHAP\_integration/R/A6H\_7\_shap\_downstream\_optimal\_params\_smll\_clusters.html
  better figures somewhere else?



**Figure 4 HiChIP**

Calling differential loops (4A)

* HiChIP/hicdc\_plus/
* filter\_loops\_for\_hicdcplus.html
* get\_counts\_for\_sig\_loops.html
* hcdc\_output\_plots.html

**Figure 5 TFs in beige loops**

MED1 signal

* HiChIP/vs\_med1/

