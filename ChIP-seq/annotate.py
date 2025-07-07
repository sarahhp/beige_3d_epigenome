"""Annotate chromatin state bed files with overlapping genes
using the ensembl genome annotation 
"""

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2022-03-24"

import os

rule annotate_cH_all:
    input:
        expand("{dir}/annotations/{sample}/{sample}_gene_filt.tab",
                dir="all_test02/model",
                sample = config["condition"].keys()),
        expand("{dir}/annotations/{sample}/{sample}_go.html",
                dir="all_test02/model",
                sample = config["condition"].keys())

rule annotate_states:
    input:
        bed = "{dir}/{sample}_20_segments.bed",
        annot = config["annot_files"]
    params:
        gene_filter =  "$14 >= 100"
    output:
        genecov = "{dir}/annotations/{sample}/{sample}_gene_filt.tab",
    shell:
        """awk '$3=="gene"' {input.annot} |
             bedtools intersect -wo -a {input.bed} -b stdin |
             awk -F "\t" '{params.gene_filter}' > {output.genecov}; 
        no_genes=$( wc -l <{output.genecov} ); 
        no_pc=$(grep -c "protein_coding" {output.genecov} );
        echo {wildcards.sample}': no. genes = '$no_genes' no. pc genes = '$no_pc;"""

rule dose_go:
    input:
        genelist = "{dir}/annotations/{sample}/{sample}_gene_filt.tab",
        background = "/projects/imb-pkbphil/sp/annotations/ensembl95/ucsc_format/min100bp_Homo_sapiens.GRCh38.95_ucsc.format.list",
        molsigdb = "/projects/imb-pkbphil/sp/annotations/molsigdb/msigdb.v7.4.symbols.gmt",
        script = "chstate_go_multi.Rmd"
    params:
        wd = os.getcwd()
    output:
        go = "{dir}/annotations/{sample}/{sample}_go.html"
    script:
        "/projects/imb-pkbphil/sp/chipseq/modelling_chstates/chromHMM/chstate_go_multi.Rmd"


#rule filter_annotation:
#    input:
#        genecov = "{dir}/{sample}_gene_intersect.tab",
#    params: 
#        gene_filter = "$14 >= 100"
#    output:
#        all = "{dir}/{sample}_gene_filt.tab",
#        pc = "{dir}/{sample}_gene_pc_filt.tab"
#    shell:
#        """awk -F "\t" '{params.gene_filter}' {input.genecov} > {output.all}; 
#        awk -F "\t" '$13 ~ "protein_coding" && {params.gene_filter}' \
#                {input.genecov}  > {output.pc}; """
#
#rule get_gene_coverage:
#   input:
#       filt = expand("{{dir}}/{sample}_gene{cat}_filt.tab",
#               sample = config["modification"]["APEX2-LMNA"],
#               cat = ["","_pc"])
#   output:
#       expand("{dir}/{dataset}_gene_coverage.tab",
#            dir = COMP_OUTDIR,
#            dataset = config["dataset"])
#   shell:
#       """for x in {input.filt}; do
#             no_genes=$( wc -l <$x )
#  
#       """
#
#

def get_samples():
  names = [f[:-13] for f in os.listdir(os.getcwd()) if f.endswith("_segments.bed")]
  return names
  
#attempt to determine number of states from file name instead of adding a config variable.  
#Not working at present 
# def is_em_file(f):
#   if f.startswith("emissions"):
#     if os.path.isfile(f):
#       if f.split(sep=".")[1] == "txt":
#         return True
#       else: 
#         return False
#     else: 
#       return False
#   else:
#     return False
# 
# def find_num_states:
#   em_files = [f for f in os.listdir(os.getcwd()) if is_em_file(f)]
#   em_file = em_files.sort()[1]
#   num_states = 15#?
#   return 15

rule heat_all:
  input:
    expand("annotations/{sample}_{coord}.{ext}",
            sample = get_samples(),
            coord = os.path.basename(config["coords"]),
            ext="plot.pdf")


rule chedarise:
  input:
    mat = expand("emissions_{states}.txt",
                states = config["num_states"]),
    script = "/projects/imb-pkbphil/sp/chipseq/modelling_chstates/chromHMM/plot_chromhmm_as_chedr.R"
  output:
    expand("emissions_{states}.plot.pdf",
              states = config["num_states"]),
    "emissions.order.txt"
  shell:
    "Rscript {input.script} {input.mat}"

rule heatmap_annotation:
  input:
    bed = "{sample}_segments.bed",
    chromHMM = "/projects/imb-pkbphil/sp/pkgs/chromHMM/ChromHMM/ChromHMM.jar"
  params:
    bin_size = 1000,
    coords_dir = config["coords"],
    out_name = lambda wildcards: "annotations/{}_{}".format(wildcards.sample,
                                    os.path.basename(config["coords"]))
  output:
    expand("annotations/{{sample}}_{coord}.{ext}", 
          ext=["png","txt","svg"],
          coord = os.path.basename(config["coords"]))
  shell:
    """
    unset DISPLAY
    java -mx8G -jar {input.chromHMM} OverlapEnrichment \
    -b {params.bin_size} {input.bed} {params.coords_dir} \
    {params.out_name}
    """

rule chedrise_enrichment:  
  input:
    mat = expand("annotations/{{sample}}_{coord}.txt",
             coord = os.path.basename(config["coords"])),
    order = "emissions.order.txt",
    script = "/projects/imb-pkbphil/sp/chipseq/modelling_chstates/chromHMM/plot_enrichments_as_chedr_maps.R"
  output:
    mat = expand("annotations/{{sample}}_{coord}.plot.pdf",
             coord = os.path.basename(config["coords"])),
  shell:
    "Rscript {input.script} {input.mat}"
