Code and analyzed data files to accompany paper Hazell Pickering et al 2025 "Multimodal remodeling of epigenetic and enhancer networks shapes the transcriptional landscape of beige adipocytes"

[https://www.biorxiv.org/content/10.1101/2025.03.28.645896v1]


**Figure 1 RNA-Seq**

Differential Expression (1E)
- RNA-seq/limma_DE_donor6_endpt.html

Deconvolution using snRNA-seq from Miranda et al.
- RNA-seq/make_custom_sigmatrix_miranda2025.html
- RNA-seq/epic_plots.html


**Figure 2 ATAC-seq**
(and chromatin states)

Chromatin states (Supp 5) 
- ChIPseq/

Gene expression in DARs (2F)
- ATAC-seq/gene_expr_in_DARs.html


**Figure 3 SHAP integration**

Normalising RNA-seq reads in promoters
- RNA-seq/promoter_counts_processing.html

Python notebooks
- SHAP_integration/py_nb/

Parameter testing and Random forest comparsion (Supp S7)
- plots_hyperparameter_search.html
- random_forest_A6H_7.html 

Optimised lightGBM model with SHAP interpretations (3B-D)
- A6H_7_features_opt_params.html
- waterfall_figures $$$ better figs elewhere? (3I,J)

Clustering SHAP values (3E-H)
- SHAP_integration/R/A6H_7_shap_downstream_optimal_params_smll_clusters.html
better figures somewhere else? 


**Figure 4 HiChIP**

Calling differential loops (4A)
- HiChIP/hicdc_plus/
 - filter_loops_for_hicdcplus.html
 - get_counts_for_sig_loops.html
 - hcdc_output_plots.html 

**Figure 5 TFs in beige loops**

MED1 signal
- HiChIP/vs_med1/

ATAC footprinting 
- to be added