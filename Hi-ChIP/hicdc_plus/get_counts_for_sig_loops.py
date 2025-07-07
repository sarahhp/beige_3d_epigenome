"""   """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2024-09-09"

config["samples"] = ["1-24770_S1","2-24771_S8","3-24772_S5","4-24773_S6"]

rule all:
    input:
        expand("hicdcplus/{sample}.sig_counts.hicdc.txt",
                sample = config["samples"]),
        "hicdcplus/all_sig_loops.txt"

rule concat_sig_loops:
    input:
        expand("{{dir}}/{sample}/{sample}.sig_loops.bedpe",
                sample = config["samples"])
    output:
         "{dir}/all_sig_loops.bed"
    shell:
         "cat {input} | sed 's/chr1\tstart1\tstart2/#chr\tstartI\tstartJ/g' | "
            "cut -f 1,2,5 | sort -n -k1.4 -k2.1 -k3.1 | uniq > {output} "

rule raw_loop_to_bed:
     input:
         "{dir}/{sample}/reg_raw.{chr}.{sample}_L001.5k.and",
         "{dir}/{sample}/reg_raw.{chr}.{sample}_L001.5k.xor"
     output:
         "{dir}/{sample}/{chr}.all_loops.bed"
     shell:
        "cat {input} | grep 'bin1_mid	bin2_mid' -v | "
            "awk -v OFS='\t' '{{ print \"{wildcards.chr}\",int($2),int($3),$4 }}' | "
            "sed '1i #chr\tstartI\tstartJ\tcount' "
            " > {output} "
        

rule merge_all_loops:
    input:
        expand("{dir}/{{sample}}/{chr}.all_loops.bed",
            dir = "maps_output",
            chr = [f'chr{i}' for i in range(1, 23)] + ["chrX"])
    output:
        "{dir}/{sample}/all_loops.bed"
    shell:
        "cat {input} | sort -n -k1.4 -k2.1 -k3.1 | uniq > {output}"
        
rule get_counts_for_sig_loops:
    input:
        counts = "maps_output/{sample}/all_loops.bed",
        sig_loops = "maps_output/all_sig_loops.bed"
    output:
        "hicdcplus/{sample}.sig_counts.hicdc.txt"
    shell:
        "bedtools intersect -loj -f 1 -r -a {input.sig_loops} -b {input.counts} | "
            "awk -v OFS='\t' '{{gsub(\"\\\.\",\"0\", $7)}} {{print $1,$2,$2+5000,$1,$3,$3+5000,$3-$2,$7}}' | "
            "sed '1i chrI\tstartI\tendI\tchrJ\tstartJ\tendJ\tD\tcounts' "
            "> {output} "

rule bed_to_text_header:
    input:
        "maps_output/all_sig_loops.bed"
    output:
        "hicdcplus/all_sig_loops.txt"
    shell:
        "sed 's/#chr/chr/g' {input} > {output} "
